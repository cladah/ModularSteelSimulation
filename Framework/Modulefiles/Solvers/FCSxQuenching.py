"""FEniCSx 0.10 thermomechanical quenching solver.

This is a cleaned-up version of the supplied code. It still depends on the
project-specific Framework.Datastream_file and FCSxHelper modules, and on the
Datastream.xdmf file being available in the working directory.
"""

from mpi4py import MPI
import numpy as np
import ufl
import basix
import dolfinx
from dolfinx import fem, mesh, default_scalar_type
from dolfinx.fem.petsc import NonlinearProblem
from scipy.spatial import cKDTree
import meshio

from Framework.Datastream_file import readdatastream, adjustdatastream
from Framework.Modulefiles.Solvers.FCSxHelper import load_datastream_to_fenicsx


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def _interp_points(V):
    """Return interpolation points across DOLFINx/Basix API variants."""
    points = V.element.interpolation_points
    return points() if callable(points) else points


def load_field_to_function(domain, filename, field_name, step_idx=0):
    """Load scalar point data from an XDMF time series into a P1 Function.

    The XDMF field must be stored as point data and its points must describe
    the same physical mesh as ``domain``. Coordinates are matched by nearest
    neighbor; this is not a substitute for validating that the meshes match.
    """
    V = fem.functionspace(domain, ("Lagrange", 1))
    function = fem.Function(V)
    function.name = field_name

    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, _cells = reader.read_points_cells()
            t, point_data, _cell_data = reader.read_data(step_idx)

        if field_name not in point_data:
            raise KeyError(
                f"Field {field_name!r} not found in {filename!r}. "
                f"Available point fields: {list(point_data.keys())}"
            )

        raw_data = np.asarray(point_data[field_name])
        if raw_data.ndim != 1:
            if raw_data.shape[-1] != 1:
                raise ValueError(
                    f"Expected scalar point data for {field_name!r}; "
                    f"got shape {raw_data.shape}."
                )
            raw_data = raw_data.reshape(-1)

        xdmf_coords = np.asarray(points)[:, : domain.geometry.dim]
        if raw_data.shape[0] != xdmf_coords.shape[0]:
            raise ValueError(
                f"Field has {raw_data.shape[0]} values but XDMF has "
                f"{xdmf_coords.shape[0]} points."
            )
    except Exception as exc:
        # Do not silently continue with a zero material field.
        raise RuntimeError(
            f"Could not load field {field_name!r} from {filename!r}"
        ) from exc

    dof_coords = V.tabulate_dof_coordinates()[:, : domain.geometry.dim]
    tree = cKDTree(xdmf_coords)
    distances, indices = tree.query(dof_coords)
    if np.any(~np.isfinite(distances)):
        raise RuntimeError(f"Invalid coordinate mapping for field {field_name!r}.")

    # Check that nearest-neighbor mapping is genuinely close to the source mesh.
    scale = max(float(np.ptp(xdmf_coords, axis=0).max()), 1.0)
    if np.max(distances, initial=0.0) > 1.0e-8 * scale:
        raise ValueError(
            f"Mesh coordinates for {field_name!r} do not match the XDMF points; "
            f"maximum nearest-neighbor distance is {np.max(distances):.3e}."
        )

    if function.x.array.size != indices.size:
        raise RuntimeError(
            "Unexpected scalar P1 DOF layout: function array and coordinate "
            "mapping have different sizes."
        )
    function.x.array[:] = np.real(raw_data[indices])
    function.x.scatter_forward()

    if domain.comm.rank == 0:
        print(f">>> Loaded {field_name!r} from step {step_idx} (t={t:.6g})")
    return function


def tensor_to_voigt(tensor, dim):
    """Convert a stress tensor to [xx, yy, xy] or [xx, yy, zz, yz, xz, xy]."""
    if dim == 2:
        return ufl.as_vector([tensor[0, 0], tensor[1, 1], tensor[0, 1]])
    if dim == 3:
        return ufl.as_vector(
            [tensor[0, 0], tensor[1, 1], tensor[2, 2],
             tensor[1, 2], tensor[0, 2], tensor[0, 1]]
        )
    raise ValueError(f"Only 2D and 3D are supported; got dim={dim}.")


def voigt_to_tensor(v, dim=2):
    """Convert [xx, yy, xy] or [xx, yy, zz, yz, xz, xy] to a tensor."""
    if dim == 2:
        return ufl.as_tensor([[v[0], v[2]], [v[2], v[1]]])
    if dim == 3:
        return ufl.as_tensor(
            [[v[0], v[5], v[4]],
             [v[5], v[1], v[3]],
             [v[4], v[3], v[2]]]
        )
    raise ValueError(f"Only 2D and 3D are supported; got dim={dim}.")


def eps(u):
    return ufl.sym(ufl.grad(u))


def sigma_TRIP(u, T, fM, minput, gdim):
    """Linear isotropic stress with thermal and martensitic eigenstrains."""
    epsilon = ufl.sym(ufl.grad(u))
    identity = ufl.Identity(gdim)

    # Model parameters retained from the supplied code.
    alpha = float(minput.get("alpha", 1.2e-5))
    T0 = float(minput.get("T0", 1113.15))
    beta = float(minput.get("beta", 0.04))
    E = float(minput.get("E", 210.0e9))
    nu = float(minput.get("nu", 0.3))

    epsilon_th = alpha * (T - T0) * identity
    epsilon_tr = (beta / 3.0) * fM * identity
    epsilon_el = epsilon - epsilon_th - epsilon_tr

    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    return lmbda * ufl.tr(epsilon_el) * identity + 2.0 * mu * epsilon_el


def TemperatureForm(domain, minput, ginput, V_T, funcs):
    """Build the transient heat residual and entropy-history objects."""
    Cv = float(minput.get("Cv", 550.0))
    Tref = float(minput.get("Tref", 293.15))
    k = float(minput.get("thermal_conductivity", 44.5))
    rho = float(minput.get("density", 7800.0))
    Tsurf = float(minput.get("Tsurf", 293.15))
    h = float(minput.get("heat_transfer_coefficient", 1000.0))
    T_initial = float(minput.get("T0", 1113.15))

    T = funcs["T"]
    dT = funcs["dT"]
    dtime = funcs["dtime"]

    # The supplied model uses s = Cv*T/Tref as its thermal history variable.
    s = (Cv / Tref) * T
    s_expr = fem.Expression(s, _interp_points(V_T))
    s_old = fem.Function(V_T)
    s_old.x.array[:] = (T_initial * Cv) / Tref
    s_old.x.scatter_forward()

    geometry_type = ginput["Geometry"]["Type"]
    if geometry_type == "Cylinder":
        radius = float(ginput["Geometry"].get("radius", 0.008))

        def boundary_marker(x):
            return np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), radius)

    elif geometry_type in ("4PointBend", "3PointBend"):
        height = float(ginput["Geometry"]["height"])

        def boundary_marker(x):
            return np.logical_or(np.isclose(x[1], height), np.isclose(x[1], 0.0))

    else:
        raise KeyError(
            f"Geometry {geometry_type!r} is not implemented in TemperatureForm."
        )

    facet_dim = domain.topology.dim - 1
    surface_facets = mesh.locate_entities_boundary(domain, facet_dim, boundary_marker)
    surface_facets = np.sort(surface_facets.astype(np.int32))
    facet_values = np.full(surface_facets.size, 1, dtype=np.int32)
    facet_tag = mesh.meshtags(domain, facet_dim, surface_facets, facet_values)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag)

    therm_res = (
        rho * Tref * (s - s_old) / dtime * dT
        + ufl.dot(k * ufl.grad(T), ufl.grad(dT))
    ) * ufl.dx + h * (T - Tsurf) * dT * ds(1)

    return therm_res, s_old, s_expr


def DisplacementForm(domain, minput, V_u, funcs):
    """Build the quasistatic mechanical residual."""
    del V_u  # The space is represented by the mixed displacement field.
    u = funcs["u"]
    du = funcs["du"]
    T = funcs["T"]
    fM = funcs["fM"]
    stress = sigma_TRIP(u, T, fM, minput, domain.geometry.dim)
    return ufl.inner(stress, eps(du)) * ufl.dx


def MartensiteForm(domain, minput, V_fM, funcs):
    """Build the algebraic martensite-fraction residual."""
    fM = funcs["fM"]
    dfM = funcs["dfM"]
    T = funcs["T"]

    Ms = load_field_to_function(domain, "Datastream.xdmf", "KM_Ms_Martensite")
    beta = load_field_to_function(domain, "Datastream.xdmf", "KM_b_Martensite")
    fMeq = float(minput.get("fMeq", 1.0))

    fM_target = ufl.conditional(
        ufl.gt(Ms, T),
        fMeq * (1.0 - ufl.exp(-beta * (Ms - T))),
        0.0,
    )
    fM_res = (fM - fM_target) * dfM * ufl.dx
    return fM_res


def _function_values_at_points(function, target_points, gdim):
    """Return Function values in target point order using nearest coordinates.

    The function must be scalar or blocked-vector P1. This is used only for
    exporting nodal fields; it does not transfer fields for solving.
    """
    V = function.function_space
    dof_coords = V.tabulate_dof_coordinates()[:, :gdim]
    block_size = V.dofmap.index_map_bs
    values = function.x.array

    if block_size > 1:
        if values.size % block_size:
            raise RuntimeError("Unexpected blocked-vector Function array size.")
        values = values.reshape((-1, block_size))
    else:
        values = values.reshape((-1,))

    if values.shape[0] != dof_coords.shape[0]:
        raise RuntimeError(
            "Cannot map function values to coordinates: DOF coordinate count "
            f"({dof_coords.shape[0]}) differs from value count ({values.shape[0]})."
        )

    tree = cKDTree(dof_coords)
    distances, indices = tree.query(target_points[:, :gdim])
    if np.any(~np.isfinite(distances)):
        raise RuntimeError("Failed to map finite-element DOFs to output points.")
    return values[indices]


# ---------------------------------------------------------------------------
# Main solver entry point
# ---------------------------------------------------------------------------
def FCSx4PB_Quench(parent):
    """Run the 4-point-bend quenching calculation and export nodal results."""
    if dolfinx.__version__.split(".")[0:2] != ["0", "10"]:
        print(f"Warning: code was targeted at DOLFINx 0.10; found {dolfinx.__version__}.")

    if MPI.COMM_WORLD.rank == 0:
        print("Using FEniCSx solver: 4-Point Bend (Quenching)")

    minput = parent.minput
    ginput = parent.ginput
    nsteps = int(minput["quench_steps"])
    if nsteps <= 0:
        raise ValueError("minput['quench_steps'] must be a positive integer.")

    domain = load_datastream_to_fenicsx("Datastream.xdmf")
    gdim = domain.geometry.dim
    if gdim not in (2, 3):
        raise ValueError(f"Only 2D and 3D meshes are supported; got gdim={gdim}.")

    xdmf_nodes = np.asarray(readdatastream("nodes"))
    if xdmf_nodes.ndim != 2 or xdmf_nodes.shape[1] < gdim:
        raise ValueError("readdatastream('nodes') did not return an (N, gdim) array.")
    xdmf_nodes = xdmf_nodes[:, :gdim]

    # Output spaces. Displacement is solved with P2 but exported as P1 nodal data.
    cell = domain.basix_cell()
    P2_u = basix.ufl.element("Lagrange", cell, 2, shape=(gdim,))
    P1_scalar = basix.ufl.element("Lagrange", cell, 1)
    P1_vector = basix.ufl.element("Lagrange", cell, 1, shape=(gdim,))

    V_el = basix.ufl.mixed_element([P2_u, P1_scalar, P1_scalar])
    V = fem.functionspace(domain, V_el)
    U = fem.Function(V)
    U.name = "solution"

    u, T, fM = ufl.split(U)
    v = ufl.TestFunction(V)
    du, dT, dfM = ufl.split(v)
    dU = ufl.TrialFunction(V)

    V_u, _ = V.sub(0).collapse()
    V_T, _ = V.sub(1).collapse()
    V_fM, _ = V.sub(2).collapse()
    V_ux, _ = V.sub(0).sub(0).collapse()
    V_uy, _ = V.sub(0).sub(1).collapse()
    V_uz = V.sub(0).sub(2).collapse()[0] if gdim == 3 else None

    dt_value = float(minput["quenchtime"]) / nsteps
    if dt_value <= 0:
        raise ValueError("The quench time must be positive.")
    dtime = fem.Constant(domain, default_scalar_type(dt_value))

    u_old = fem.Function(V_u)
    T_old = fem.Function(V_T)
    fM_old = fem.Function(V_fM)
    T_initial = float(minput.get("T0", 1113.15))
    T_old.x.array[:] = T_initial
    fM_old.x.array[:] = 0.0
    T_old.x.scatter_forward()
    fM_old.x.scatter_forward()

    funcs = {
        "u": u, "du": du, "u_old": u_old,
        "T": T, "dT": dT, "T_old": T_old,
        "fM": fM, "dfM": dfM, "fM_old": fM_old,
        "dtime": dtime,
    }

    # Boundary conditions: retain the constraints from the supplied model.
    def b1_BC(x):
        return np.logical_and(np.isclose(x[0], 0.0), np.isclose(x[1], 0.0))

    def sym_x_BC(x):
        return np.isclose(x[0], 0.0)
    def sym_thickness(x):
        return np.isclose(x[2], 0.0)

    ux_value = fem.Function(V_ux)
    uy_value = fem.Function(V_uy)
    uz_value = fem.Function(V_uz)
    ux_value.x.array[:] = 0.0
    uy_value.x.array[:] = 0.0
    uz_value.x.array[:] = 0.0

    dofs_y = fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), b1_BC)
    dofs_x = fem.locate_dofs_geometrical((V.sub(0).sub(0), V_ux), sym_x_BC)
    dofs_z = fem.locate_dofs_geometrical((V.sub(0).sub(2), V_uz), sym_thickness)
    bcs = [
        fem.dirichletbc(uy_value, dofs_y, V.sub(0).sub(1)),
        fem.dirichletbc(ux_value, dofs_x, V.sub(0).sub(0)),
        fem.dirichletbc(uz_value, dofs_z, V.sub(0).sub(2)),
    ]

    if gdim == 3:
        def sym_z_BC(x):
            return np.isclose(x[2], 0.0)

        uz_value = fem.Function(V_uz)
        uz_value.x.array[:] = 0.0
        dofs_z = fem.locate_dofs_geometrical((V.sub(0).sub(2), V_uz), sym_z_BC)
        bcs.append(fem.dirichletbc(uz_value, dofs_z, V.sub(0).sub(2)))

    fM_res = MartensiteForm(domain, minput, V_fM, funcs)
    T_res, s_old, s_expr = TemperatureForm(domain, minput, ginput, V_T, funcs)
    u_res = DisplacementForm(domain, minput, V_u, funcs)

    residual = u_res + T_res + fM_res
    jacobian = ufl.derivative(residual, U, dU)

    petsc_options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",
        "snes_atol": 1.0e-7,
        "snes_rtol": 1.0e-7,
        "snes_max_it": 50,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }
    problem = NonlinearProblem(
        residual,
        U,
        bcs=bcs,
        J=jacobian,
        petsc_options=petsc_options,
        petsc_options_prefix="martensite_solver_",
    )

    if MPI.COMM_WORLD.rank == 0:
        print(">>> Starting solve...")

    for step in range(nsteps):
        if MPI.COMM_WORLD.rank == 0:
            print(f"Step {step + 1}/{nsteps}")

        problem.solve()
        U.x.scatter_forward()
        reason = problem.solver.getConvergedReason()
        if reason <= 0:
            raise RuntimeError(
                f"Nonlinear solver failed at step {step + 1}; PETSc SNES reason={reason}."
            )

        # Update history fields from the converged mixed solution.
        u_old.interpolate(U.sub(0))
        T_old.interpolate(U.sub(1))
        fM_old.interpolate(U.sub(2))
        s_old.interpolate(s_expr)
        u_old.x.scatter_forward()
        T_old.x.scatter_forward()
        fM_old.x.scatter_forward()
        s_old.x.scatter_forward()

    # Collapse mixed fields and interpolate displacement from P2 to P1 for output.
    u_final = U.sub(0).collapse()
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    V_out_scalar = fem.functionspace(domain, ("Lagrange", 1))
    V_out_vector = fem.functionspace(domain, ("Lagrange", 1, (gdim,)))
    V_out_stress = fem.functionspace(
        domain, ("Lagrange", 1, (3 if gdim == 2 else 6,))
    )

    u_out = fem.Function(V_out_vector, name="Displacement")
    u_out.interpolate(u_final)
    T_out = fem.Function(V_out_scalar, name="Temperature")
    T_out.interpolate(T_final)
    fM_out = fem.Function(V_out_scalar, name="Martensite")
    fM_out.interpolate(fM_final)

    stress_ufl = sigma_TRIP(U.sub(0), U.sub(1), U.sub(2), minput, gdim)
    stress_voigt_ufl = tensor_to_voigt(stress_ufl, dim=gdim)
    stress_expr = fem.Expression(stress_voigt_ufl, _interp_points(V_out_stress))
    stress_out = fem.Function(V_out_stress, name="Stress")
    stress_out.interpolate(stress_expr)
    stress_out.x.scatter_forward()

    # Map each output function by its own DOF coordinates, not by assuming
    # finite-element DOFs have the same ordering as geometry nodes.
    displacement_values = _function_values_at_points(u_out, xdmf_nodes, gdim)
    temperature_values = _function_values_at_points(T_out, xdmf_nodes, gdim)
    martensite_values = _function_values_at_points(fM_out, xdmf_nodes, gdim)
    stress_values = _function_values_at_points(stress_out, xdmf_nodes, gdim)

    result_dict = {
        "Displacement": displacement_values,
        "Temperature": temperature_values,
        "Martensite": martensite_values,
        "Austenite": 1.0 - martensite_values,
        "Stress": np.real(stress_values),
    }
    adjustdatastream(result_dict, datapos="nodes", t_data=float(minput["quenchtime"]))

    if MPI.COMM_WORLD.rank == 0:
        print(">>> Quenching solve completed and results exported.")

    return result_dict
