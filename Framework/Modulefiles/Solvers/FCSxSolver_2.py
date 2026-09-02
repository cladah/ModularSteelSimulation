from mpi4py import MPI
from dolfinx import mesh, fem, io, nls
import ufl
import basix
from Framework.Datastream_file import readdatastream, readDatastreamMesh, adjustdatastream
from dolfinx.fem.petsc import NonlinearProblem
from scipy.interpolate import LinearNDInterpolator
import pyvista as pv
from dolfinx.io import XDMFFile
from dolfinx import fem, mesh, plot, default_scalar_type
from basix.ufl import element, mixed_element
from petsc4py import PETSc
from scipy.spatial import cKDTree
import meshio
import dolfinx
import numpy as np


def tensor_to_voigt(tensor, dim=2):
    if dim == 2:
        # Standard 2D Voigt (Plane Stress/Strain)
        return ufl.as_vector([tensor[0, 0],
                              tensor[1, 1],
                              tensor[0, 1]])
    else:
        # Standard 3D Voigt
        return ufl.as_vector([tensor[0, 0],
                              tensor[1, 1],
                              tensor[2, 2],
                              tensor[1, 2],
                              tensor[0, 2],
                              tensor[0, 1]])

def voigt_to_tensor(v, dim=2):
    if dim == 2:
        # Vector: [sigma_xx, sigma_yy, sigma_xy]
        return ufl.as_tensor([[v[0], v[2]],
                              [v[2], v[1]]])
    else:
        # Vector: [sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy]
        return ufl.as_tensor([[v[0], v[5], v[4]],
                              [v[5], v[1], v[3]],
                              [v[4], v[3], v[2]]])

def load_datastream_to_fenicsx(filename="Datastream.xdmf"):
    """
    Reads the mesh from Datastream.xdmf and converts it
    into a DOLFINx (FEniCSx) Mesh/Domain.
    """
    try:
        # 1. Open the file in read mode
        # We use the COMM_WORLD communicator for parallel compatibility
        with XDMFFile(MPI.COMM_WORLD, filename, "r") as xdmf:

            # 2. Read the mesh (domain)
            # 'mesh' is the default name meshio uses in the XDMF header
            domain = xdmf.read_mesh(name="mesh")

        print(f"Successfully loaded {filename} into FEniCSx domain.")
        print(
            f"Nodes: {domain.geometry.x.shape[0]}, Cells: {domain.topology.index_map(domain.topology.dim).size_global}")

        return domain

    except Exception as e:
        print(f"Error loading mesh to FEniCSx: {e}")
        return None

def load_field_to_function(domain, filename, field_name, step_idx=0):
    """
    Loads a specific data field from an XDMF using meshio.TimeSeriesReader
    and injects it into a FEniCSx Function.
    """
    # 1. Create the FunctionSpace and Function
    # Most XDMF node data is stored as Lagrange Degree 1 (linear)
    V = dolfinx.fem.functionspace(domain, ("Lagrange", 1))
    u = dolfinx.fem.Function(V)
    u.name = field_name


    # 2. Use TimeSeriesReader to get the raw NumPy data
    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            # Read metadata for the requested step
            t, point_data, cell_data = reader.read_data(step_idx)
            if field_name not in point_data:
                raise KeyError(f"Field '{field_name}' not found in {filename}. "
                               f"Available: {list(point_data.keys())}")

            raw_data = point_data[field_name]
            raw_data = raw_data.flatten()
            xdmf_coords = points[:, :domain.geometry.dim]
    except Exception as e:
        print(f">>> Error reading {filename} with meshio: {e}")
        return u  # Returns zero-initialized function as fallback

    dof_coords = V.tabulate_dof_coordinates()[:, :domain.geometry.dim]
    tree = cKDTree(xdmf_coords)
    _, xdmf_indices = tree.query(dof_coords)

    # 3. Inject data into the FEniCSx Function
    # NOTE: This assumes the mesh in 'domain' has the same node ordering as 'filename'
    # .real handles potential complex-type builds of PETSc
    if u.x.array.shape == raw_data.shape:
        u.x.array[:] = raw_data[xdmf_indices].real
    else:
        # If there is a mismatch (e.g. ghost cells or MPI), we only fill the local nodes
        # This part is crucial for stability in parallel/larger meshes
        local_size = V.dofmap.index_map.size_local
        u.x.array[:local_size] = raw_data[xdmf_indices[:local_size]].real

    # Update ghost values across MPI processes
    u.x.scatter_forward()

    print(f">>> Successfully loaded '{field_name}' from step {step_idx} (t={t:.2e}s)")
    return u

def voigt_strain(symm_tensor):
    return ufl.as_vector([symm_tensor[0, 0], symm_tensor[1, 1], 2.0 * symm_tensor[0, 1]])

def eps(u):
    return ufl.sym(ufl.grad(u))

def sigma(u, T, fM, minput, gdim):
    E, nu = 210e9, 0.3
    if minput["Materialtype"] == "Elastic":
        return sigma_elastic()
    elif minput["Materialtype"] == "LinearPlastic":
        return sigma_TRIP()

def sigma_elastic(elastic_eps, minput):
    E, nu = 210e9, 0.3
    mu = E / (2 * (1 + nu))
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lmbda * ufl.tr(elastic_eps) * ufl.Identity(len(elastic_eps)) + 2.0 * mu * elastic_eps

def sigma_TRIP(u, T, fM, minput, gdim):
    # 1. Kinematics
    epsilon = ufl.sym(ufl.grad(u))
    I = ufl.Identity(gdim)

    # 2. Thermal Strain: epsilon_th = alpha * (T - T0) * I
    # (Assuming alpha and T0 are in minput)
    alpha = 1.2e-5
    T0 = 1113.15
    epsilon_th = alpha * (T - T0) * I

    # 3. Transformation Strain: 3% volumetric expansion -> 1% linear
    # epsilon_tr = (beta / 3) * fM * I
    beta = 0.04
    epsilon_tr = (beta / 3.0) * fM * I

    # 4. Elastic Strain
    epsilon_el = epsilon - epsilon_th - epsilon_tr

    # 5. Hooke's Law (Lame Parameters)
    E, nu = 210e9, 0.3
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    lmbda = E * nu / (1 - nu**2)
    mu = E / (2 * (1 + nu))

    return lmbda * ufl.tr(epsilon_el) * I + 2 * mu * epsilon_el
def vM(s):
    epsilon = 1e-6
    sb = s - (1. / 3.) * ufl.tr(s) * ufl.Identity(ufl.shape(s)[0])
    return ufl.sqrt(3. / 2. * ufl.inner(sb, sb) + epsilon)


def lode_func(sigma):
    """
    Computes the Lode angle (theta) from a stress tensor using UFL.
    Returns a value where cos(3*theta) is bounded between -1 and 1.
    """
    # 1. Identity tensor and Deviatoric stress
    d = ufl.shape(sigma)[0]
    I = ufl.Identity(d)
    s = sigma - (1 / d) * ufl.tr(sigma) * I

    # 2. Invariants
    J2 = 0.5 * ufl.inner(s, s)
    # J3 = det(s)
    J3 = ufl.det(s)

    # 3. Prevent division by zero at hydrostatic states
    eps = 1e-12

    # 4. Calculate cos(3*theta)
    cos3theta = (3 * ufl.sqrt(3) / 2) * (J3 / (J2 ** 1.5 + eps))

    # 5. Clamp the value to [-1, 1] to avoid NaNs in acos due to precision
    cos3theta = ufl.conditional(ufl.gt(cos3theta, 1.0), 1.0, cos3theta)
    cos3theta = ufl.conditional(ufl.lt(cos3theta, -1.0), -1.0, cos3theta)
    #L = -cos(3 * theta)
    # Return theta
    return cos3theta#(1 / 3) * ufl.acos(cos3theta)
def triax(sig):
    return sig_h(sig)/vM(sig)

def sigma_von_mises(u, minput):
    """Calculates the von Mises stress expression using UFL."""
    epsilon = 1e-6
    s = sigma(u, minput)
    sb = s - (1. / 3.) * ufl.tr(s) * ufl.Identity(len(u))
    return ufl.sqrt(3. / 2. * ufl.inner(sb, sb) + epsilon)

def sig_h(sigma):
    return (1. / 3.) * ufl.tr(sigma)

def MartensiteForm(domain, minput, V_fM, funcs):
    fM, fM_old, dfM = funcs["fM"], funcs["fM_old"], funcs["dfM"]
    T, u, dtime = funcs["T"], funcs["u"], funcs["dtime"]

    # Load parameters from the datastream we created/interpolated earlier
    Ms = load_field_to_function(domain, "Datastream.xdmf", "KM_Ms_Martensite")
    beta = load_field_to_function(domain, "Datastream.xdmf", "KM_b_Martensite")
    fMeq = 1.0
    Ms_expr = Ms + 1e-7 * sigma_von_mises(u,minput)
    # Koistinen-Marburger Expression for current state
    fM_target = ufl.conditional(ufl.gt(Ms_expr, T),
                                fMeq * (1.0 - ufl.exp(-beta * (Ms_expr - T))),
                                0.0)

    fM_expr = fem.Expression(fM_target, V_fM.element.interpolation_points)

    # Residual: Standard L2 Projection of the algebraic KM equation
    fM_res = (fM - fM_target) * dfM * ufl.dx

    return fM_res, fM_expr

def DisplacementForm(domain, minput, V_u, funcs):
    wtC_data = readdatastream("Composition_C")
    fM = funcs["fM"]
    fM_old = funcs["fM_old"]
    dfM = funcs["dfM"]
    u = funcs["u"]
    du = funcs["du"]
    T = funcs["T"]
    T_old = funcs["T_old"]
    dT = funcs["dT"]
    #s = sigma(u, minput)
    gdim = domain.geometry.dim
    s = sigma_TRIP(u, T, fM, minput, gdim)
    disp_res = ufl.inner(s, eps(du)) * ufl.dx
    return disp_res


def DisplacementForm_temp(domain, minput, V_u, funcs):
    # 1. Extract necessary functions
    u = funcs["u"]  # Current displacement
    du = funcs["du"]  # Test function
    T = funcs["T"]  # Current temperature
    fM = funcs["fM"]

    # 2. Define Thermal Parameters
    # alpha is the thermal expansion coefficient (e.g., 1.2e-5 for steel)
    # T_ref is the stress-free reference temperature
    alpha = minput.get("alpha", 1.2e-5)
    T_ref = 293.15
    dim = domain.topology.dim

    # 3. Define Strains
    def eps(v):
        return ufl.sym(ufl.grad(v))

    # Thermal strain: alpha * deltaT * Identity
    # ufl.Identity(dim) creates a 2x2 or 3x3 identity matrix
    eps_th = alpha * (T - T_ref) * ufl.Identity(dim)
    eps_fM = 0.01 * fM * ufl.Identity(dim)
    # 4. Define Stress (Constitutive Equation)
    # Use (Total Strain - Thermal Strain) to get Elastic Strain
    def sigma_thermal(v, temp, fM):
        elastic_strain = eps(v) - eps_th + eps_fM
        # Use your existing sigma logic, but pass the elastic_strain
        # If your sigma() helper is already defined, you may need to
        # modify it to accept the strain tensor directly.
        return sigma_elastic(elastic_strain, minput)

    # 5. Residual Formulation
    # Inner product of stress and the variation of strain
    disp_res = ufl.inner(sigma_thermal(u, T, fM), eps(du)) * ufl.dx

    return disp_res

def LargeDispForm(domain, minput, V_u, funcs):
    u = funcs["u"]
    u_old = funcs["u_old"]
    du = funcs["du"]

    # Spatial dimension
    d = len(u)

    # Identity tensor
    I = ufl.variable(ufl.Identity(d))

    # Deformation gradient
    F = ufl.variable(I + ufl.grad(u))

    # Right Cauchy-Green tensor
    C = ufl.variable(F.T * F)

    # Invariants of deformation tensors
    Ic = ufl.variable(ufl.tr(C))
    J = ufl.variable(ufl.det(F))

    E = (1.0e4)
    nu = (0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))

    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2

def TemperatureForm(domain, minput, ginput, V_T, funcs):
    Cv = 550.0 - 8.314
    T0 = 1113.15
    Tref = 293.15
    T = funcs["T"]
    dT = funcs["dT"]
    u = funcs["u"]
    dtime = funcs["dtime"]
    k = 44.5
    rho = 7800
    Tsurf = 293.15

    s = Cv / Tref * T
    s_expr = fem.Expression(s, V_T.element.interpolation_points)
    s_old = fem.Function(V_T)
    s_old.x.array[:] = (T0 * Cv) / Tref
    if ginput["Geometry"]["Type"] == "Cylinder":
        f_map = lambda x: np.isclose(np.sqrt(x[0]**2+x[1]**2), 0.008)
    elif ginput["Geometry"]["Type"] in ["4PointBend","3PointBend"]:
        f_map = lambda x: np.logical_or(np.isclose(x[1], ginput["Geometry"]["height"]),np.isclose(x[1], 0.0))
    else:
        raise KeyError("Geometry not implemented in Temperature form")

    surface_facets = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, f_map)
    facet_tag = mesh.meshtags(domain, domain.topology.dim - 1, surface_facets, np.full(len(surface_facets), 1, dtype=np.int32))
    surf_ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag)
    therm_res = (rho * Tref * (s - s_old) / dtime * dT + ufl.dot(k * ufl.grad(T), ufl.grad(dT))) * ufl.dx \
                + 1000.0 * (T - Tsurf) * dT * surf_ds(1)

    return therm_res, s_old, s_expr

def plotMesh(domain):
    topology, cell_types, geometry = plot.vtk_mesh(domain, domain.topology.dim)
    grid = pv.UnstructuredGrid(topology, cell_types, geometry)

    plotter = pv.Plotter()
    plotter.add_mesh(grid, show_edges=True, edge_color="black", color="lightblue")

    plotter.add_axes()
    plotter.set_background("white")
    print(">>> Opening mesh preview window...")
    plotter.show(title="Imported Mesh Preview")

def FCSx4PB_Force(parent):
    print('Using FEniCSx solver: 4-Point Bend (Displacement/Force Control)')
    minput = parent.minput
    ginput = parent.ginput

    # 1. Load Mesh
    domain = load_datastream_to_fenicsx("Datastream.xdmf")
    gdim = domain.geometry.dim
    #plotMesh(domain)

    # 2. Map meshio nodes to FEniCSx nodes for result extraction
    xdmf_nodes = readdatastream("nodes")  # (N, 2) or (N, 3)
    fenics_nodes = domain.geometry.x[:, :gdim]
    from scipy.spatial import cKDTree
    tree = cKDTree(fenics_nodes)
    _, fenics_to_xdmf_map = tree.query(xdmf_nodes[:, :gdim])

    # 3. Setup Spaces
    P2_vect = basix.ufl.element("Lagrange", domain.basix_cell(), 2, shape=(gdim,))
    P1 = basix.ufl.element("Lagrange", domain.basix_cell(), 1)

    voigt_dim = 3 if gdim == 2 else 6
    V_P1 = fem.functionspace(domain, ("Lagrange", 1, (1,)))
    V_P1_voigt = fem.functionspace(domain, ("Lagrange", 1, (voigt_dim,)))
    V_P1_vct = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    V_P2 = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    V_P2_voight = fem.functionspace(domain, ("Lagrange", 2, (voigt_dim,)))
    V_P2_vct = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))

    # 2. Combine into a Mixed Element
    V_el = basix.ufl.mixed_element([P2_vect, P1, P1])
    V = fem.functionspace(domain, V_el)

    # Functions for solving
    U = fem.Function(V)
    u, T, fM = ufl.split(U)
    V_test = ufl.TestFunction(V)
    du, dT, dfM = ufl.split(V_test)
    dU = ufl.TrialFunction(V)

    # Subspaces for BCs and interpolation
    V_u, _ = V.sub(0).collapse()
    V_ux, _ = V.sub(0).sub(0).collapse()
    V_uy, _ = V.sub(0).sub(1).collapse()
    V_uz, _ = V.sub(0).sub(2).collapse()
    V_T, _ = V.sub(1).collapse()
    V_fM, _ = V.sub(2).collapse()

    # 4. History and Constants
    u_old = fem.Function(V_u)
    #T_old = fem.Function(V_T)
    #T_old.x.array[:] = 1113.15
    #T_old.x.array[:] = 273.15 + 20
    #fM_old = fem.Function(V_fM)
    T_old = load_field_to_function(domain, "Datastream.xdmf", "Temperature")
    fM_old = load_field_to_function(domain, "Datastream.xdmf", "Martensite")

    funcs = {"u": u, "du": du, "u_old": u_old, "T": T, "dT": dT,
             "T_old": T_old, "fM": fM, "dfM": dfM, "fM_old": fM_old}

    # 5. Boundary Conditions (4-Point Bending)
    def b1_BC(x): return np.logical_and(np.isclose(x[0], 0.0), np.isclose(x[1], 0.0))

    if ginput["Geometry"]["Type"] == "4PointBend":
        def t1_BC(x): return np.logical_and(np.isclose(x[0], ginput["Geometry"]["width"]/3), np.isclose(x[1], ginput["Geometry"]["height"]))

    elif ginput["Geometry"]["Type"] == "3PointBend":
        def t1_BC(x):
            return np.logical_and(np.isclose(x[0], ginput["Geometry"]["width"]/2), np.isclose(x[1], ginput["Geometry"]["height"]))
    else:
        raise KeyError("Geometry type not implemented in Bending solver")

    def sym_x_BC(x): return np.isclose(x[0], ginput["Geometry"]["width"]/2)
    def sym_z_BC(x): return np.isclose(x[2], 0.00)

    # Dirichlet setup
    uP_val = fem.Function(V_uy)
    uP_val.x.array[:] = 0.0
    uzero_val = fem.Function(V_ux)
    uzero_val.x.array[:] = 0.0
    u_zero_func = fem.Function(V_u)
    u_zero_func.x.array[:] = 0.0

    bcs = [
        fem.dirichletbc(uP_val, fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), t1_BC), V.sub(0).sub(1)),
        fem.dirichletbc(uzero_val, fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), b1_BC), V.sub(0).sub(1)),
        fem.dirichletbc(uzero_val, fem.locate_dofs_geometrical((V.sub(0).sub(0), V_ux), sym_x_BC), V.sub(0).sub(0)),
        fem.dirichletbc(uzero_val, fem.locate_dofs_geometrical((V.sub(0).sub(2), V_uz), sym_z_BC), V.sub(0).sub(2))
    ]

    # 6. Variational Forms
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    #T_res, s_old, s_expr = TemperatureForm(domain, minput, V_T, funcs)
    T_dummy = (T - T_old) * dT * ufl.dx
    fM_dummy = (fM - fM_old) * dfM * ufl.dx
    u_res = DisplacementForm(domain, minput, V_u, funcs)

    #Res = u_res + T_dummy + fM_dummy
    Res = u_res + fM_res + T_dummy
    Jac = ufl.derivative(Res, U, dU)

    # 7. Solver Setup
    uP_val.x.array[:] = -minput["Umax"]

    petsc_options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "none",
        "snes_atol": 1e-4,
        "snes_rtol": 1e-4,
        "snes_monitor": None,
        "ksp_error_if_not_converged": True,
        "ksp_type": "gmres",
        "ksp_rtol": 1e-7,
        "ksp_monitor": None,
        "pc_type": "hypre",
        "pc_hypre_type": "boomeramg",
        "pc_hypre_boomeramg_max_iter": 1,
        "pc_hypre_boomeramg_cycle_type": "v",
    }
    problem = NonlinearProblem(Res, U, bcs=bcs, J=Jac, petsc_options=petsc_options, petsc_options_prefix="Bending_")

    print(">>> Starting Solve...")
    for step in range(minput["loadsteps"]):
        print(f"Step nr {step + 1}")
        problem.solve()
        converged = problem.solver.getConvergedReason()
        num_iter = problem.solver.getIterationNumber()
        assert converged > 0, f"Solver did not converge, got {converged}."
        print(
            f"Solver converged after {num_iter} iterations with converged reason {converged}."
        )
        U.x.scatter_forward()

        # Update History
        T_old.interpolate(U.sub(1))
        u_old.interpolate(U.sub(0))
        #s_old.interpolate(s_expr)

    u_final = U.sub(0).collapse()
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    def project_and_map(ufl_expr, target_space):
        expr = fem.Expression(ufl_expr, target_space.element.interpolation_points)
        f = fem.Function(target_space)
        f.interpolate(expr)
        return f

    vm_stress = project_and_map(sigma_von_mises(U.sub(0), minput), V_P1)
    sh = project_and_map(sig_h(sigma(U.sub(0), minput)), V_P1)
    trax = project_and_map(triax(sigma(U.sub(0), minput)), V_P1)
    lode = project_and_map(lode_func(sigma(U.sub(0), minput)), V_P1)

    stress_ufl = tensor_to_voigt(sigma(U.sub(0), minput), dim=gdim)
    stress_out = project_and_map(stress_ufl, V_P1_voigt)

    u_out = fem.Function(V_P1_vct)
    u_out.interpolate(U.sub(0))

    eps_ufl = tensor_to_voigt(eps(U.sub(0)), dim=gdim)
    eps_out = project_and_map(eps_ufl, V_P1_voigt)

    print(np.shape(vm_stress.x.array))
    print(vm_stress.x.array)
    print(np.max(vm_stress.x.array))
    print(np.min(vm_stress.x.array))

    res_dict = {
        "Displacement": u_out.x.array.reshape(-1, gdim)[fenics_to_xdmf_map],
        "Strain_el": eps_out.x.array.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Temperature": T_final.x.array[fenics_to_xdmf_map],
        "vonMises": vm_stress.x.array[fenics_to_xdmf_map],
        "Martensite": fM_final.x.array[fenics_to_xdmf_map],
        "Austenite": 1.0 - fM_final.x.array[fenics_to_xdmf_map],
        "Stress": stress_out.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Stress_hydrostatic": sh.x.array[fenics_to_xdmf_map],
        "Triaxiality": trax.x.array[fenics_to_xdmf_map],
        "Lode": lode.x.array[fenics_to_xdmf_map]
    }

    adjustdatastream(res_dict, datapos="nodes", t_data=0.0)
    #plot_4PB_results(U, minput)
def FCSx4PB_Quench(parent):
    print('Using FEniCSx solver: 4-Point Bend (Quenching)')
    minput = parent.minput
    ginput = parent.ginput
    T0 = 1113.15

    # 1. Load Mesh
    domain = load_datastream_to_fenicsx("Datastream.xdmf")
    gdim = domain.geometry.dim
    # plotMesh(domain)

    # 2. Map meshio nodes to FEniCSx nodes for result extraction
    xdmf_nodes = readdatastream("nodes")  # (N, 2) or (N, 3)
    fenics_nodes = domain.geometry.x[:, :gdim]
    from scipy.spatial import cKDTree
    tree = cKDTree(fenics_nodes)
    _, fenics_to_xdmf_map = tree.query(xdmf_nodes[:, :gdim])

    # 3. Setup Spaces
    P2 = basix.ufl.element("Lagrange", domain.basix_cell(), 2, shape=(domain.geometry.dim,))
    P1 = basix.ufl.element("Lagrange", domain.basix_cell(), 1)

    voigt_dim = 3 if gdim == 2 else 6
    V_P1 = fem.functionspace(domain, ("Lagrange", 1, (1,)))
    V_P1_voigt = fem.functionspace(domain, ("Lagrange", 1, (voigt_dim,)))
    V_P1_vct = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    V_P2 = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    V_P2_vct = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))

    # 2. Combine into a Mixed Element
    V_el = basix.ufl.mixed_element([P2, P1, P1])

    # 3. Create the FunctionSpace using the lowercase factory function
    V = fem.functionspace(domain, V_el)

    # Functions for solving
    U = fem.Function(V)
    u, T, fM = ufl.split(U)
    V_test = ufl.TestFunction(V)
    du, dT, dfM = ufl.split(V_test)
    dU = ufl.TrialFunction(V)


    # Subspaces for BCs and interpolation
    V_u, _ = V.sub(0).collapse()
    V_ux, _ = V.sub(0).sub(0).collapse()
    V_uy, _ = V.sub(0).sub(1).collapse()
    V_uz, _ = V.sub(0).sub(2).collapse()
    V_T, _ = V.sub(1).collapse()
    V_fM, _ = V.sub(2).collapse()

    # 4. History and Constants
    dtime = fem.Constant(domain, default_scalar_type(minput["quenchtime"] / minput["quench_steps"]))
    u_old = fem.Function(V_u)
    T_old = fem.Function(V_T)
    T_old.x.array[:] = 1113.15
    fM_old = fem.Function(V_fM)






    funcs = {"u": u, "du": du, "u_old": u_old, "T": T, "dT": dT,
             "T_old": T_old, "fM": fM, "dfM": dfM, "fM_old": fM_old, "dtime": dtime}

    # 5. Boundary Conditions (4-Point Bending)
    def b1_BC(x):
        return np.logical_and(np.isclose(x[0], 0.0), np.isclose(x[1], 0.0))

    def sym_x_BC(x):
        return np.isclose(x[0], ginput["Geometry"]["width"]/2)

    def sym_z_BC(x):
        return np.isclose(x[2], 0.00)

    # Dirichlet setup
    ux_sym_val = fem.Function(V_ux)
    ux_sym_val.x.array[:] = 0.0
    uz_sym_val = fem.Function(V_uz)
    uz_sym_val.x.array[:] = 0.0
    u_zero_func = fem.Function(V_u)
    u_zero_func.x.array[:] = 0.0

    bcs = [
        fem.dirichletbc(ux_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), b1_BC), V.sub(0).sub(1)),
        fem.dirichletbc(ux_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(0), V_ux), sym_x_BC), V.sub(0).sub(0)),
        fem.dirichletbc(uz_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(2), V_uz), sym_z_BC), V.sub(0).sub(2)),
    ]

    # 6. Variational Forms
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    T_res, s_old, s_expr = TemperatureForm(domain, minput, ginput, V_T, funcs)
    T_dummy = (T - T_old) * dT * ufl.dx
    fM_dummy = (fM - fM_old) * dfM * ufl.dx
    u_res = DisplacementForm(domain, minput, V_u, funcs)

    Res = u_res + T_res + fM_res
    Jac = ufl.derivative(Res, U, dU)
    petsc_options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",  # Backtracking: prevents overshooting
        "snes_atol": 1e-7,  # Absolute tolerance
        "snes_rtol": 1e-7,  # Relative tolerance
        "snes_max_it": 50,  # Allow more iterations for difficult steps
        "ksp_type": "preonly",  # Direct solver doesn't need KSP iterations
        "pc_type": "lu",  # LU decomposition
        "pc_factor_mat_solver_type": "mumps",  # High-performance direct solver
    }
    # 7. Solver Setup
    problem = NonlinearProblem(Res, U, bcs=bcs, J=Jac, petsc_options=petsc_options, petsc_options_prefix="martensite_solver_")

    # 8. Time Loop
    print(">>> Starting Solve...")
    for step in range(minput["quench_steps"]):
        print(f"Step nr {step + 1}")

        problem.solve()
        U.x.scatter_forward()

        # Update History
        T_old.interpolate(U.sub(1))
        u_old.interpolate(U.sub(0))
        s_old.interpolate(s_expr)
        # fM_old.interpolate(fem.Expression(fM_exp, V_fM.element.interpolation_points()))
    # 9. Extract Results for Datastream
    # Collapse to extract arrays
    u_final = U.sub(0).collapse()
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    vm_expr_ufl = sigma_von_mises(U.sub(0), minput)
    vm_expr = fem.Expression(vm_expr_ufl, V_P1.element.interpolation_points)
    vm_stress = fem.Function(V_P1)
    vm_stress.interpolate(vm_expr)

    lode_expr_ufl = lode_func(sigma(U.sub(0), minput))
    lode_expr = fem.Expression(lode_expr_ufl, V_P1.element.interpolation_points)
    lode = fem.Function(V_P1)
    lode.interpolate(lode_expr)

    trax_expr_ufl = triax(sigma(U.sub(0), minput))
    trax_expr = fem.Expression(trax_expr_ufl, V_P1.element.interpolation_points)
    trax = fem.Function(V_P1)
    trax.interpolate(trax_expr)

    sh_expr_ufl = sig_h(sigma(U.sub(0), minput))
    sh_expr = fem.Expression(sh_expr_ufl, V_P1.element.interpolation_points)
    sh = fem.Function(V_P1)
    sh.interpolate(sh_expr)

    s_expr_ufl = sigma(U.sub(0), minput)
    sv_expr_ufl = tensor_to_voigt(s_expr_ufl, dim=gdim)
    stress_expr = fem.Expression(sv_expr_ufl, V_P1_voigt.element.interpolation_points)
    stress = fem.Function(V_P1_voigt)
    stress.interpolate(stress_expr)

    # Use the map to export in the correct meshio order

    res_dict = {
        "Displacement": u_final.x.array.reshape(-1, gdim)[fenics_to_xdmf_map],
        "Temperature": T_final.x.array[fenics_to_xdmf_map],
        "vonMises": vm_stress.x.array[fenics_to_xdmf_map],
        "Martensite": fM_final.x.array[fenics_to_xdmf_map],
        "Austenite": 1.0 - fM_final.x.array[fenics_to_xdmf_map],
        "Stress": stress.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Stress_hydrostatic": sh.x.array[fenics_to_xdmf_map],
        "Triaxiality": trax.x.array[fenics_to_xdmf_map],
        "Lode": lode.x.array[fenics_to_xdmf_map]
    }

    adjustdatastream(res_dict, datapos="nodes", t_data=0.0)

def Cylinder_2D_Quench(parent):
    print('Using FEniCSx solver: 2D Cylinder (Quenching)')
    minput = parent.minput
    ginput = parent.ginput
    T0 =  1113.15

    # 1. Load Mesh
    domain = load_datastream_to_fenicsx("Datastream.xdmf")
    gdim = domain.geometry.dim

    # 2. Map meshio nodes to FEniCSx nodes for result extraction
    xdmf_nodes = readdatastream("nodes")  # (N, 2) or (N, 3)
    xdmf_nodes_2d = xdmf_nodes[:, :gdim]
    fenics_nodes = domain.geometry.x[:, :gdim]
    from scipy.spatial import cKDTree
    tree = cKDTree(fenics_nodes)
    distance, fenics_to_xdmf_map = tree.query(xdmf_nodes_2d)

    #_, fenics_to_xdmf_map = tree.query(xdmf_nodes[:, :gdim])

    # 3. Setup Spaces
    P2 = basix.ufl.element("Lagrange", domain.basix_cell(), 2, shape=(domain.geometry.dim,))
    P1 = basix.ufl.element("Lagrange", domain.basix_cell(), 1)

    voigt_dim = 3 if gdim == 2 else 6
    V_P1 = fem.functionspace(domain, ("Lagrange", 1, (1,)))
    V_P1_voigt = fem.functionspace(domain, ("Lagrange", 1, (voigt_dim,)))
    V_P1_vct = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    V_P2 = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    V_P2_vct = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))

    # 2. Combine into a Mixed Element
    V_el = basix.ufl.mixed_element([P2, P1, P1])

    # 3. Create the FunctionSpace using the lowercase factory function
    V = fem.functionspace(domain, V_el)

    # Functions for solving
    U = fem.Function(V)
    u, T, fM = ufl.split(U)
    V_test = ufl.TestFunction(V)
    du, dT, dfM = ufl.split(V_test)
    dU = ufl.TrialFunction(V)

    # Subspaces for BCs and interpolation
    V_u, _ = V.sub(0).collapse()
    V_ux, _ = V.sub(0).sub(0).collapse()
    V_uy, _ = V.sub(0).sub(1).collapse()
    V_T, _ = V.sub(1).collapse()
    V_fM, _ = V.sub(2).collapse()

    # 4. History and Constants
    dtime = fem.Constant(domain, default_scalar_type(minput["quenchtime"] / minput["quench_steps"]))
    u_old = fem.Function(V_u)
    T_old = fem.Function(V_T)
    T_old.x.array[:] = T0
    fM_old = fem.Function(V_fM)

    funcs = {"u": u, "du": du, "u_old": u_old, "T": T, "dT": dT,
             "T_old": T_old, "fM": fM, "dfM": dfM, "fM_old": fM_old, "dtime": dtime}

    # 5. Boundary Conditions (Cylinder)

    def sym_x_BC(x): return np.isclose(x[0], 0.0)

    def Surface_T_BC(x): return np.isclose(np.sqrt(x[0]**2 + x[1]**2), 0.00)

    def sym_y_BC(x): return np.isclose(x[1], 0.00)

    # Dirichlet setup
    uT_val = fem.Function(V_uy)
    uT_val.x.array[:] = T0
    ux_sym_val = fem.Function(V_ux)
    ux_sym_val.x.array[:] = 0.0
    uy_sym_val = fem.Function(V_uy)
    uy_sym_val.x.array[:] = 0.0
    u_zero_func = fem.Function(V_u)
    u_zero_func.x.array[:] = 0.0

    bcs = [
        #fem.dirichletbc(uT_val, fem.locate_dofs_geometrical((V.sub(1), V_T), Surface_T_BC), V.sub(1)),
        fem.dirichletbc(ux_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(0), V_ux), sym_x_BC), V.sub(0).sub(0)),
        fem.dirichletbc(uy_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), sym_y_BC), V.sub(0).sub(1))
    ]

    # 6. Variational Forms
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    T_res, s_old, s_expr = TemperatureForm(domain, minput, ginput, V_T, funcs)
    u_res = DisplacementForm(domain, minput, V_u, funcs)
    T_dummy = (T - T_old) * dT * ufl.dx
    fM_dummy = (fM - fM_old) * dfM * ufl.dx
    #u_dummy = (u - u_old) * du * ufl.dx

    Res = u_res + T_res + fM_res
    Jac = ufl.derivative(Res, U, dU)

    # 1. Define the problem with a unique prefix
    use_superlu = PETSc.IntType == np.int64  # or PETSc.ScalarType == np.complex64
    sys = PETSc.Sys()  # type: ignore
    if sys.hasExternalPackage("mumps") and not use_superlu:
        linear_solver = "mumps"
    elif sys.hasExternalPackage("superlu_dist"):
        linear_solver = "superlu_dist"
    else:
        linear_solver = "petsc"
    petsc_options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "none",
        "snes_atol": 1e-4,
        "snes_rtol": 1e-4,
        "snes_monitor": None,
        "ksp_error_if_not_converged": True,
        "ksp_type": "gmres",
        "ksp_rtol": 1e-7,
        "ksp_monitor": None,
        "pc_type": "hypre",
        "pc_hypre_type": "boomeramg",
        "pc_hypre_boomeramg_max_iter": 1,
        "pc_hypre_boomeramg_cycle_type": "v",
    }
    petsc_options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",  # Backtracking: prevents overshooting
        "snes_atol": 1e-7,  # Absolute tolerance
        "snes_rtol": 1e-7,  # Relative tolerance
        "snes_max_it": 50,  # Allow more iterations for difficult steps
        "ksp_type": "preonly",  # Direct solver doesn't need KSP iterations
        "pc_type": "lu",  # LU decomposition
        "pc_factor_mat_solver_type": "mumps",  # High-performance direct solver
    }

    # 7. Solver Setup
    problem = NonlinearProblem(Res, U, bcs=bcs, J=Jac, petsc_options=petsc_options, petsc_options_prefix="Cylinder_")
    #problem = NonlinearProblem(Res, U, bcs=bcs, J=Jac , petsc_options_prefix="martensite_solver_")
    # 8. Time Loop
    print(">>> Starting Solve...")
    for step in range(minput["quench_steps"]):
        print(f"Step nr {step + 1}")
        #uT_val.x.array[:] = 273.15 + 20

        problem.solve()
        U.x.scatter_forward()

        # Update History
        T_old.interpolate(U.sub(1))
        u_old.interpolate(U.sub(0))
        s_old.interpolate(s_expr)

    # 9. Extract Results for Datastream
    # Collapse to extract arrays
    u_final = U.sub(0).collapse()
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    vm_expr_ufl = sigma_von_mises(U.sub(0), minput)
    vm_expr = fem.Expression(vm_expr_ufl, V_P1.element.interpolation_points)
    vm_stress = fem.Function(V_P1)
    vm_stress.interpolate(vm_expr)

    lode_expr_ufl = lode_func(sigma(U.sub(0),minput))
    lode_expr = fem.Expression(lode_expr_ufl, V_P1.element.interpolation_points)
    lode = fem.Function(V_P1)
    lode.interpolate(lode_expr)

    trax_expr_ufl = triax(sigma(U.sub(0), minput))
    trax_expr = fem.Expression(trax_expr_ufl, V_P1.element.interpolation_points)
    trax = fem.Function(V_P1)
    trax.interpolate(trax_expr)

    sh_expr_ufl = sig_h(sigma(U.sub(0), minput))
    sh_expr = fem.Expression(sh_expr_ufl, V_P1.element.interpolation_points)
    sh = fem.Function(V_P1)
    sh.interpolate(sh_expr)

    s_expr_ufl = sigma(U.sub(0), minput)
    sv_expr_ufl = tensor_to_voigt(s_expr_ufl, dim=gdim)
    stress_expr = fem.Expression(sv_expr_ufl, V_P1_voigt.element.interpolation_points)
    stress = fem.Function(V_P1_voigt)
    stress.interpolate(stress_expr)

    s_expr_ufl = sigma_TRIP(U.sub(0), U.sub(1), U.sub(2), minput, gdim)
    sv_expr_ufl = tensor_to_voigt(s_expr_ufl, dim=gdim)
    stress_expr = fem.Expression(sv_expr_ufl, V_P1_voigt.element.interpolation_points)
    stress = fem.Function(V_P1_voigt)
    stress.interpolate(stress_expr)

    eps_total = ufl.sym(ufl.grad(U.sub(0)))

    # 2. Convert to Voigt Notation (to store as a vector: [exx, eyy, exy] for 2D)
    ev_expr_ufl = tensor_to_voigt(eps_total, dim=gdim)
    strain_expr = fem.Expression(ev_expr_ufl, V_P1_voigt.element.interpolation_points)
    strain_fn = fem.Function(V_P1_voigt)
    strain_fn.interpolate(strain_expr)

    alpha = minput.get("alpha", 1.2e-5)
    beta = 0.03
    I = ufl.Identity(gdim)
    eps_th = alpha * (U.sub(1) - T0) * I
    eps_tr = (beta / 3.0) * U.sub(2) * I

    # 3. Elastic Strain Tensor
    eps_el_ufl = eps_total - eps_th - eps_tr

    # Interpolate Elastic Strain (Voigt)
    ev_el_expr_ufl = tensor_to_voigt(eps_el_ufl, dim=gdim)
    eps_el_expr = fem.Expression(ev_el_expr_ufl, V_P1_voigt.element.interpolation_points)
    eps_el_fn = fem.Function(V_P1_voigt)
    eps_el_fn.interpolate(eps_el_expr)

    # Interpolate Elastic Strain (Voigt)
    ev_tr_expr_ufl = tensor_to_voigt(eps_tr, dim=gdim)
    eps_tr_expr = fem.Expression(ev_tr_expr_ufl, V_P1_voigt.element.interpolation_points)
    eps_tr_fn = fem.Function(V_P1_voigt)
    eps_tr_fn.interpolate(eps_tr_expr)

    # Interpolate Elastic Strain (Voigt)
    ev_th_expr_ufl = tensor_to_voigt(eps_th, dim=gdim)
    eps_th_expr = fem.Expression(ev_th_expr_ufl, V_P1_voigt.element.interpolation_points)
    eps_th_fn = fem.Function(V_P1_voigt)
    eps_th_fn.interpolate(eps_th_expr)

    res_dict = {
        "Displacement": u_final.x.array.reshape(-1, gdim)[fenics_to_xdmf_map],
        "Temperature": T_final.x.array[fenics_to_xdmf_map],
        "vonMises": vm_stress.x.array[fenics_to_xdmf_map],
        "Martensite": fM_final.x.array[fenics_to_xdmf_map],
        "Austenite": 1.0 - fM_final.x.array[fenics_to_xdmf_map],
        "Stress": stress.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Strain": strain_fn.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Strain_el": eps_el_fn.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Strain_th": eps_th_fn.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Strain_trip": eps_tr_fn.x.array.real.reshape(-1, voigt_dim)[fenics_to_xdmf_map],
        "Stress_hydrostatic": sh.x.array[fenics_to_xdmf_map],
        "Triaxiality": trax.x.array[fenics_to_xdmf_map],
        "Lode": lode.x.array[fenics_to_xdmf_map]
    }

    adjustdatastream(res_dict, datapos="nodes", t_data=0.0)

def plot_4PB_results(U,minput):
    """Helper to visualize the Temperature and Martensite fields."""
    import pyvista as pv

    domain = U.function_space.mesh
    gdim = domain.geometry.dim
    voigt_dim = 3 if gdim == 2 else 6

    # 1. Extract sub-functions
    u_P2 = U.sub(0).collapse()
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    # 2. Create a Linear Plotting Space for Displacement (P1)
    V_P1 = fem.functionspace(domain, ("Lagrange", 1, (1,)))
    V_P1_vec = fem.functionspace(domain, ("Lagrange", 1, (gdim,)))
    V_P1_voigt = fem.functionspace(domain, ("Lagrange", 1, (voigt_dim,)))

    u_P1 = fem.Function(V_P1_vec)
    u_P1.interpolate(u_P2)  # Map P2 data down to P1 nodes

    # 3. Create PyVista Grid
    topology, cell_types, geometry = plot.vtk_mesh(domain, gdim)
    grid = pv.UnstructuredGrid(topology, cell_types, geometry)

    # 4. Calculate Displacement Magnitude from the P1 function
    u_values = u_P1.x.array.real.reshape(-1, gdim)
    u_magnitude = np.linalg.norm(u_values, axis=1)

    s_expr_ufl = sigma(U.sub(0), minput)
    sv_expr_ufl = tensor_to_voigt(s_expr_ufl, dim=gdim)
    stress_expr = fem.Expression(sv_expr_ufl, V_P1_voigt.element.interpolation_points())
    stress = fem.Function(V_P1_voigt)
    stress.interpolate(stress_expr)

    vm_expr_ufl = sigma_von_mises(U.sub(0), minput)
    vm_expr = fem.Expression(vm_expr_ufl, V_P1.element.interpolation_points())
    vm_stress = fem.Function(V_P1)
    vm_stress.interpolate(vm_expr)

    # 5. Attach data to grid
    grid.point_data["Displacement"] = u_values
    grid.point_data["Displacement Magnitude"] = u_magnitude
    grid.point_data["Temperature"] = T_final.x.array.real
    grid.point_data["Martensite"] = fM_final.x.array.real
    grid.point_data["vonMises"] = vm_stress.x.array
    warped_grid = grid.warp_by_vector("Displacement", factor=1)

    # 6. Set up Plotter
    p = pv.Plotter(shape=(1, 2), window_size=[1200, 400])

    p.subplot(0, 0)
    p.add_text("Displacement Magnitude", font_size=10)
    p.add_mesh(warped_grid, scalars="Displacement Magnitude", cmap="jet", show_edges=False)
    p.view_xy()
    p.show_grid()

    p.subplot(0, 1)
    p.add_text("von-Mises stress", font_size=10)
    p.add_mesh(grid, scalars="vonMises", cmap="jet", show_edges=False, clim=[0.0, 1.0e9])
    p.view_xy()
    p.show_grid(font_size=10)
    p.show()

    p2 = pv.Plotter(shape=(1, 2), window_size=[1200, 400])

    p2.subplot(0, 0)
    p2.add_text("Displacement Magnitude", font_size=10)
    p2.add_mesh(warped_grid, scalars="Displacement Magnitude", cmap="jet", show_edges=False)
    p2.view_xy()
    p2.show_grid()

    p2.subplot(0, 1)
    p2.add_text("Martensite Fraction", font_size=10)
    p2.add_mesh(grid, scalars="Martensite", cmap="turbo", show_edges=False, clim=[0.0, 1.0])
    p2.view_xy()

    """
    p.subplot(0, 2)
    p.add_text("Final Temperature", font_size=10)
    p.add_mesh(grid, scalars="Temperature", cmap="inferno", show_edges=True)
    p.view_xy()
    """
    p2.show()
def plot_quenching_results(U):
    """Helper to visualize the Temperature and Martensite fields."""
    import pyvista as pv

    # Extract sub-functions
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    # Create PyVista Grid
    topology, cell_types, geometry = plot.vtk_mesh(T_final.function_space)
    grid = pv.UnstructuredGrid(topology, cell_types, geometry)

    # Attach data
    grid.point_data["Temperature"] = T_final.x.array.real
    grid.point_data["Martensite"] = fM_final.x.array.real

    # Set up Plotter
    p = pv.Plotter(shape=(1, 2))

    p.subplot(0, 0)
    p.add_text("Final Temperature")
    p.add_mesh(grid, scalars="Temperature", cmap="inferno", show_edges=True)

    p.subplot(0, 1)
    p.add_text("Martensite Fraction")
    p.add_mesh(grid, scalars="Martensite", cmap="viridis", show_edges=True)

    p.show()

