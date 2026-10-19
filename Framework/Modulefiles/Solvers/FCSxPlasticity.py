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
from ..Solvers.FCSxHelper import load_field_to_function, tensor_to_voigt, voigt_to_tensor, load_datastream_to_fenicsx, \
    ufl_piecewise_linear


def MartensiteForm(domain, minput, V_fM, funcs):
    fM, fM_old, dfM = funcs["fM"], funcs["fM_old"], funcs["dfM"]
    T, u = funcs["T"], funcs["u"]

    # Load parameters from the datastream we created/interpolated earlier
    Ms = load_field_to_function(domain, "Datastream.xdmf", "KM_Ms_Martensite")
    beta = load_field_to_function(domain, "Datastream.xdmf", "KM_b_Martensite")
    fMeq = 1.0

    fM_target = ufl.conditional(ufl.gt(Ms, T),
                                fMeq * (1.0 - ufl.exp(-beta * (Ms - T))),
                                0.0)
    fM_expr = fem.Expression(fM_target, V_fM.element.interpolation_points)

    # Residual: Standard L2 Projection of the algebraic KM equation
    fM_res = (fM - fM_target) * dfM * ufl.dx

    return fM_res, fM_expr

def TemperatureForm(domain, minput, ginput, V_T, funcs):
    Cv = 550.0 - 8.314
    T0 = minput.get("CNtemp", 1113.15)
    Tref = minput.get("Tref", 293.15)
    T = funcs["T"]
    dT = funcs["dT"]
    u = funcs["u"]
    dtime = funcs["dtime"]
    k = 30.5
    rho = 7800
    Tsurf = 293.15

    s = Cv / Tref * T
    s_expr = fem.Expression(s, V_T.element.interpolation_points)
    s_old = fem.Function(V_T)
    s_old.x.array[:] = (T0 * Cv) / Tref
    if ginput["Geometry"]["Type"] == "Cylinder":
        radius = ginput.get("radius", 0.008)
        f_map = lambda x: np.isclose(np.sqrt(x[0]**2+x[1]**2), radius)
    elif ginput["Geometry"]["Type"] in ["4PointBend", "3PointBend"]:
        f_map = lambda x: np.logical_or(np.isclose(x[1], ginput["Geometry"]["height"]),np.isclose(x[1], 0.0))
    else:
        raise KeyError("Geometry not implemented in Temperature form")

    surface_facets = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, f_map)
    facet_tag = mesh.meshtags(domain, domain.topology.dim - 1, surface_facets, np.full(len(surface_facets), 1, dtype=np.int32))
    surf_ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag)
    therm_res = (rho * Tref * (s - s_old) / dtime * dT + ufl.dot(k * ufl.grad(T), ufl.grad(dT))) * ufl.dx \
                + 1000.0 * (T - Tsurf) * dT * surf_ds(1)

    return therm_res, s_old, s_expr

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


def Cylinder_2D_Quench(parent):
    """
    2D cylinder quenching with:

        - FEniCSx nonlinear thermo-mechanical solve
        - Martensitic transformation
        - Small-strain J2 plasticity
        - Isotropic linear hardening
        - Quadrature-point plastic history variables
        - Radial-return mapping

    IMPORTANT:
        The 2D mesh is treated as PLANE STRAIN.

        Therefore the constitutive model is internally 3D:
            xx, yy, zz, xy, yz, xz

        with:
            eps_zz = 0
            eps_yz = 0
            eps_xz = 0

        Thermal and transformation strains are applied isotropic

    Plastic history variables:
        eps_p_old : quadrature-point plastic strain tensor
        ep_eq_old : quadrature-point equivalent plastic strain
    """

    print(
        "Using FEniCSx solver: "
        "2D Cylinder (Quenching + Quadrature J2 Plasticity)"
    )

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    minput = parent.minput
    ginput = parent.ginput

    T0 = minput.get("CNtemp",1113.15)

    domain = load_datastream_to_fenicsx("Datastream.xdmf")
    assert domain.comm.size == size

    gdim = domain.geometry.dim
    # ============================================================
    # 2. XDMF -> FEniCSx node mapping
    # ============================================================

    xdmf_nodes = readdatastream("nodes")

    xdmf_nodes_2d = xdmf_nodes[:, :gdim]

    fenics_nodes = (domain.geometry.x[:, :gdim])

    from scipy.spatial import cKDTree

    tree = cKDTree(fenics_nodes)

    distance, fenics_to_xdmf_map = tree.query(
        xdmf_nodes_2d
    )

    if np.max(distance) > 1e-8:
        print(
            "WARNING: XDMF/FEniCS node mapping "
            f"maximum distance = {np.max(distance):.3e}"
        )

    # ============================================================
    # 3. Function spaces
    # ============================================================

    P2 = basix.ufl.element(
        "Lagrange",
        domain.basix_cell(),
        2,
        shape=(gdim,)
    )

    P1 = basix.ufl.element(
        "Lagrange",
        domain.basix_cell(),
        1
    )

    V_P1 = fem.functionspace(
        domain,
        ("Lagrange", 1, (1,))
    )

    V_P1_vct = fem.functionspace(
        domain,
        ("Lagrange", 1, (gdim,))
    )

    voigt_dim_2d = 3

    V_P1_voigt = fem.functionspace(
        domain,
        ("Lagrange", 1, (voigt_dim_2d,))
    )

    V_P2_vct = fem.functionspace(
        domain,
        ("Lagrange", 2, (gdim,))
    )

    # ============================================================
    # 4. Mixed displacement / temperature / martensite space
    # ============================================================

    V_el = basix.ufl.mixed_element([P2, P1, P1])

    V = fem.functionspace(domain, V_el)

    U = fem.Function(V)

    u, T, fM = ufl.split(U)

    V_test = ufl.TestFunction(V)

    du, dT, dfM = ufl.split(V_test)

    dU = ufl.TrialFunction(V)

    # ============================================================
    # 5. Collapsed spaces
    # ============================================================

    V_u, _ = V.sub(0).collapse()

    V_ux, _ = (V.sub(0).sub(0).collapse())

    V_uy, _ = (V.sub(0).sub(1).collapse())
    V_uz, _ = (V.sub(0).sub(2).collapse())

    V_T, _ = V.sub(1).collapse()

    V_fM, _ = V.sub(2).collapse()

    # ============================================================
    # 6. Time step
    # ============================================================

    dtime = fem.Constant(domain,
        default_scalar_type(minput["quenchtime"] / minput["quench_steps"]))

    # ============================================================
    # 7. Previous global fields
    # ============================================================

    u_old = fem.Function(V_u)

    T_old = fem.Function(V_T)

    T_old.x.array[:] = T0

    fM_old = fem.Function(V_fM)

    # ============================================================
    # 8. Material parameters
    # ============================================================

    E = minput.get(
        "E",
        210.0e9
    )

    nu = minput.get(
        "nu",
        0.30
    )

    sigma_y0 = minput.get(
        "yield_stress",
        minput.get(
            "sigma_y",
            1200.0e6
        )
    )

    H = minput.get(
        "hardening_modulus",
        minput.get(
            "H",
            1.0e9
        )
    )

    alpha = minput.get(
        "alpha",
        1.2e-5
    )

    beta = minput.get(
        "beta",
        0.03
    )

    # ============================================================
    # 9. Lame constants
    # ============================================================

    mu = E / (
        2.0 * (1.0 + nu)
    )

    lam = (
        E * nu
        / (
            (1.0 + nu)
            * (1.0 - 2.0 * nu)
        )
    )

    # ============================================================
    # 10. Quadrature degree
    # ============================================================

    q_degree = int(
        minput.get(
            "plasticity_quadrature_degree",
            2
        )
    )

    # ============================================================
    # 11. Quadrature elements
    #
    # Plastic strain has 6 components:
    #
    # [xx, yy, zz, xy, yz, xz]
    #
    # Equivalent plastic strain is scalar.
    # ============================================================

    q_cell = domain.basix_cell()

    Q_scalar_el = basix.ufl.quadrature_element(
        q_cell,
        degree=q_degree
    )

    Q_tensor_el = basix.ufl.quadrature_element(
        q_cell,
        degree=q_degree,
        value_shape=(6,)
    )

    Q_scalar = fem.functionspace(
        domain,
        Q_scalar_el
    )

    Q_tensor = fem.functionspace(
        domain,
        Q_tensor_el
    )

    # ============================================================
    # 12. Plastic history variables
    #
    # THESE ARE QUADRATURE-POINT VARIABLES.
    # ============================================================

    eps_p_old = fem.Function(
        Q_tensor
    )

    ep_eq_old = fem.Function(
        Q_scalar
    )

    eps_p_old.x.array[:] = 0.0

    ep_eq_old.x.array[:] = 0.0

    # ============================================================
    # 13. Tensor helpers
    # ============================================================

    def voigt6_to_tensor(v):
        """
        [xx, yy, zz, xy, yz, xz]
        ->
        3x3 symmetric tensor
        """

        return ufl.as_tensor([
            [v[0], v[3], v[5]],
            [v[3], v[1], v[4]],
            [v[5], v[4], v[2]]
        ])

    def tensor_to_voigt6(t):
        return ufl.as_vector([
            t[0, 0],
            t[1, 1],
            t[2, 2],
            t[0, 1],
            t[1, 2],
            t[0, 2]
        ])

    def tensor_to_voigt2(t):
        """
        In-plane output only:
            [xx, yy, xy]
        """

        return ufl.as_vector([
            t[0, 0],
            t[1, 1],
            t[0, 1]
        ])

    # ============================================================
    # 14. Identity
    # ============================================================

    I3 = ufl.Identity(3)

    # ============================================================
    # 15. Plane-strain total strain
    #
    # u is 2D:
    #
    # [ux]
    # [uy]
    #
    # Therefore:
    #
    # eps_zz = 0
    # eps_xz = 0
    # eps_yz = 0
    # ============================================================

    eps_2d = ufl.sym(
        ufl.grad(u)
    )

    eps_total = ufl.as_tensor([
        [
            eps_2d[0, 0],
            eps_2d[0, 1],
            0.0
        ],
        [
            eps_2d[1, 0],
            eps_2d[1, 1],
            0.0
        ],
        [
            0.0,
            0.0,
            0.0
        ]
    ])

    # ============================================================
    # 16. Thermal strain
    # ============================================================

    eps_th = (
        alpha
        * (T - T0)
        * I3
    )

    # ============================================================
    # 17. Transformation strain
    # ============================================================

    eps_tr = (
        beta / 3.0
    ) * fM * I3

    # ============================================================
    # 18. Previous plastic strain tensor
    # ============================================================

    eps_p_old_tensor = (
        voigt6_to_tensor(
            eps_p_old
        )
    )

    # ============================================================
    # 19. Trial elastic strain
    # ============================================================

    eps_el_trial = (
        eps_total
        - eps_th
        - eps_tr
        - eps_p_old_tensor
    )

    # ============================================================
    # 20. Elastic trial stress
    # ============================================================

    sigma_trial = (
        lam
        * ufl.tr(eps_el_trial)
        * I3
        + 2.0
        * mu
        * eps_el_trial
    )

    # ============================================================
    # 21. Trial deviatoric stress
    # ============================================================

    s_trial = (
        sigma_trial
        - (
            ufl.tr(sigma_trial)
            / 3.0
        ) * I3
    )

    # ============================================================
    # 22. Trial von-Mises stress
    # ============================================================

    q_trial = ufl.sqrt(
        1.5
        * ufl.inner(
            s_trial,
            s_trial
        )
        + 1.0e-20
    )

    # ============================================================
    # 23. Isotropic hardening
    # ============================================================

    sigma_y_old = (sigma_y0 + H * ep_eq_old)

    # ============================================================
    # 24. Trial yield function
    # ============================================================

    f_trial = (q_trial - sigma_y_old)

    # ============================================================
    # 25. Radial return
    # ============================================================

    dgamma = ufl.conditional(ufl.gt(f_trial, 0.0),
                             f_trial / (3.0 * mu + H), 0.0)

    # ============================================================
    # 26. Flow direction
    # ============================================================

    n_flow = (1.5 * s_trial / q_trial)

    # ============================================================
    # 27. Updated plastic strain
    # ============================================================

    eps_p_new_tensor = (eps_p_old_tensor + dgamma * n_flow)

    # ============================================================
    # 28. Updated equivalent plastic strain
    # ============================================================

    ep_eq_new = (ep_eq_old + dgamma)

    # ============================================================
    # 29. Updated elastic strain
    # ============================================================

    eps_el_new = (eps_total - eps_th
                  - eps_tr - eps_p_new_tensor)

    # ============================================================
    # 30. Updated stress
    # ============================================================

    sigma = (lam * ufl.tr(eps_el_new) * I3
             + 2.0 * mu * eps_el_new)

    # ============================================================
    # 31. Plasticity state dictionary
    # ============================================================

    funcs = {
        "u": u,
        "du": du,
        "u_old": u_old,

        "T": T,
        "dT": dT,
        "T_old": T_old,

        "fM": fM,
        "dfM": dfM,
        "fM_old": fM_old,

        "dtime": dtime,

        # Plasticity
        "eps_p_old": eps_p_old,
        "ep_eq_old": ep_eq_old,

        "eps_p_new": eps_p_new_tensor,
        "ep_eq_new": ep_eq_new,

        "sigma_trial": sigma_trial,
        "sigma": sigma,

        "dgamma": dgamma,
        "yield_function": f_trial
    }

    # ============================================================
    # 32. Boundary conditions
    # ============================================================

    def sym_x_BC(x):
        return np.isclose(x[0], 0.0)

    def sym_y_BC(x):
        return np.isclose(x[1], 0.0)

    ux_sym_val = fem.Function(V_ux)
    ux_sym_val.x.array[:] = 0.0

    uy_sym_val = fem.Function(V_uy)
    uy_sym_val.x.array[:] = 0.0

    uz_zero = fem.Function(V_uz)
    uz_zero.x.array[:] = 0.0

    bcs = [
        fem.dirichletbc(ux_sym_val,
                        fem.locate_dofs_geometrical((V.sub(0).sub(0), V_ux), sym_x_BC),
                        V.sub(0).sub(0)),

        fem.dirichletbc(uy_sym_val,
                        fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), sym_y_BC),
                        V.sub(0).sub(1)),

        fem.dirichletbc(uz_zero,
                        fem.locate_dofs_geometrical((V.sub(0).sub(2), V_uz),
                                                    lambda x: np.full(x.shape[1], True, dtype=bool)),
                        V.sub(0).sub(2))
    ]

    # ============================================================
    # 33. Thermal and martensitic residuals
    # ============================================================

    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)

    T_res, s_old, s_expr = TemperatureForm(domain, minput, ginput, V_T, funcs)

    # ============================================================
    # 34. Mechanical residual
    # ============================================================

    u_res = (ufl.inner(sigma, ufl.sym(ufl.grad(du))) * ufl.dx)

    # ============================================================
    # 35. Global residual
    # ============================================================

    Res = (
        u_res
        + T_res
        + fM_res
    )

    # ============================================================
    # 36. Automatic Jacobian
    # ============================================================

    Jac = ufl.derivative(Res, U, dU)

    # ============================================================
    # 37. PETSc solver
    # ============================================================
    max_newton_iterations = 10
    petsc_options = {

        "snes_type":
            "newtonls",

        "snes_linesearch_type":
            "bt",

        "snes_atol":
            1.0e-7,

        "snes_rtol":
            1.0e-7,

        "snes_max_it":
            max_newton_iterations,

        "snes_monitor":
            None,

        "ksp_type":
            "preonly",

        "pc_type":
            "lu",

        "pc_factor_mat_solver_type":
            "mumps"
    }

    problem = NonlinearProblem(
        Res,
        U,
        bcs=bcs,
        J=Jac,
        petsc_options=petsc_options,
        petsc_options_prefix="CylinderJ2_"
    )

    print(
        ">>> Starting "
        "quadrature-point J2 plasticity solve..."
    )

    t = 0.0
    t_end = float(minput["quenchtime"])
    dt = t_end / float(minput["quench_steps"])
    dt_min = float(minput.get("dt_min", dt / 10000.0))
    dt_max = float(minput.get("dt_max", dt * 1000.0))

    max_retries = int(minput.get("max_step_retries", 20))
    target_iterations = int(minput.get("target_newton_iterations", 8))
    step = 0

    # Timestep growth / reduction factors
    dt_growth = float(minput.get("dt_growth", 1.5))

    dt_cut = float(minput.get("dt_cut", 0.5))

    while t < t_end - 1.0e-14:
        dt = min(dt, t_end - t)
        dtime.value = dt

        U_backup = U.x.array.copy()

        u_old_backup = u_old.x.array.copy()
        T_old_backup = T_old.x.array.copy()
        fM_old_backup = fM_old.x.array.copy()

        eps_p_old_backup = eps_p_old.x.array.copy()
        ep_eq_old_backup = ep_eq_old.x.array.copy()

        s_old_backup = s_old.x.array.copy()

        converged = False

        for retry in range(max_retries + 1):

            # Restore state before each attempt
            U.x.array[:] = U_backup
            U.x.scatter_forward()

            u_old.x.array[:] = u_old_backup
            T_old.x.array[:] = T_old_backup
            fM_old.x.array[:] = fM_old_backup

            eps_p_old.x.array[:] = eps_p_old_backup
            ep_eq_old.x.array[:] = ep_eq_old_backup

            s_old.x.array[:] = s_old_backup

            # Update actual timestep
            dtime.value = dt

            try:
                problem.solve()

                U.x.scatter_forward()
                newton_iterations = problem.solver.getIterationNumber()

                if newton_iterations >= max_newton_iterations:
                    raise Exception("To many newton iterations")

                converged = True
            except Exception as exc:

                print("")
                print("!!! Nonlinear solve not converged.")

                print(f"Reason: {exc}")
                print("")

                converged = False

            if converged:
                print(
                    ">>> Nonlinear solve converged."
                )

                break
            dt *= dt_cut

            print(
                f">>> Reducing timestep to "
                f"{dt:.6e}"
            )
            if dt < dt_min:
                raise RuntimeError(
                    "Adaptive timestep became smaller than "
                    f"dt_min = {dt_min:.6e}."
                )

                # Make sure we do not overshoot final time
            dt = min(
                dt,
                t_end - t
            )

        if not converged:
            raise RuntimeError(
                "Adaptive timestep failed to converge."
            )

        # --------------------------------------------------------
        # Update global history
        # --------------------------------------------------------

        u_old.interpolate(U.sub(0))
        u_old.x.scatter_forward()

        T_old.interpolate(U.sub(1))
        T_old.x.scatter_forward()

        fM_old.interpolate(U.sub(2))
        fM_old.x.scatter_forward()

        # ========================================================
        # Update quadrature plastic history
        #
        # eps_p_new_tensor is a UFL expression evaluated at
        # quadrature points.
        # ========================================================

        eps_p_new_voigt = (
            tensor_to_voigt6(
                eps_p_new_tensor
            )
        )

        eps_p_update_expr = fem.Expression(
            eps_p_new_voigt,
            Q_tensor.element.interpolation_points
        )

        eps_p_update = fem.Function(
            Q_tensor
        )

        eps_p_update.interpolate(
            eps_p_update_expr
        )

        eps_p_old.x.array[:] = (
            eps_p_update.x.array
        )
        eps_p_old.x.scatter_forward()

        # ========================================================
        # Equivalent plastic strain history
        # ========================================================

        ep_eq_update_expr = fem.Expression(
            ep_eq_new,
            Q_scalar.element.interpolation_points
        )

        ep_eq_update = fem.Function(
            Q_scalar
        )

        ep_eq_update.interpolate(
            ep_eq_update_expr
        )

        ep_eq_old.x.array[:] = (
            ep_eq_update.x.array
        )
        ep_eq_old.x.scatter_forward()

        # --------------------------------------------------------
        # Update thermal history used by your existing
        # TemperatureForm
        # --------------------------------------------------------

        s_old.interpolate(s_expr)
        s_old.x.scatter_forward()

        t += dt
        step += 1

        print(
            f">>> ACCEPTED: t = {t:.6e}, "
            f"dt = {dt:.6e}"
        )
        dt = min(dt * dt_growth, dt_max)
        dt = min(dt, t_end - t)
    # ============================================================
    # 40. Final solution fields
    # ============================================================

    u_final = U.sub(0).collapse()

    T_final = U.sub(1).collapse()

    fM_final = U.sub(2).collapse()

    # ============================================================
    # 41. Reconstruct final constitutive state
    # ============================================================

    eps_2d_final = ufl.sym(
        ufl.grad(
            U.sub(0)
        )
    )

    eps_total_final = ufl.as_tensor([
        [
            eps_2d_final[0, 0],
            eps_2d_final[0, 1],
            0.0
        ],
        [
            eps_2d_final[1, 0],
            eps_2d_final[1, 1],
            0.0
        ],
        [
            0.0,
            0.0,
            0.0
        ]
    ])

    eps_th_final = (
        alpha
        * (U.sub(1) - T0)
        * I3
    )

    eps_tr_final = (
        beta / 3.0
    ) * U.sub(2) * I3

    eps_p_final_tensor = (
        voigt6_to_tensor(
            eps_p_old
        )
    )

    eps_el_final = (
        eps_total_final
        - eps_th_final
        - eps_tr_final
        - eps_p_final_tensor
    )

    sigma_final = (
        lam
        * ufl.tr(eps_el_final)
        * I3
        + 2.0
        * mu
        * eps_el_final
    )

    # ============================================================
    # 42. Final von-Mises stress
    # ============================================================

    s_final = (
        sigma_final
        - (
            ufl.tr(sigma_final)
            / 3.0
        ) * I3
    )

    vm_final = ufl.sqrt(
        1.5
        * ufl.inner(
            s_final,
            s_final
        )
        + 1.0e-20
    )

    vm_expr = fem.Expression(
        vm_final,
        V_P1.element.interpolation_points
    )

    vm_stress = fem.Function(
        V_P1
    )

    vm_stress.interpolate(
        vm_expr
    )

    # ============================================================
    # 43. Hydrostatic stress
    # ============================================================

    sh_final = (
        ufl.tr(sigma_final)
        / 3.0
    )

    sh_expr = fem.Expression(
        sh_final,
        V_P1.element.interpolation_points
    )

    sh = fem.Function(
        V_P1
    )

    sh.interpolate(
        sh_expr
    )

    # ============================================================
    # 44. Stress triaxiality
    # ============================================================

    trax_final = (
        sh_final
        / (
            vm_final
            + 1.0e-20
        )
    )

    trax_expr = fem.Expression(
        trax_final,
        V_P1.element.interpolation_points
    )

    trax = fem.Function(
        V_P1
    )

    trax.interpolate(
        trax_expr
    )

    # ============================================================
    # 45. Lode parameter
    # ============================================================

    lode_expr_ufl = lode_func(
        sigma_final
    )

    lode_expr = fem.Expression(
        lode_expr_ufl,
        V_P1.element.interpolation_points
    )

    lode = fem.Function(
        V_P1
    )

    lode.interpolate(
        lode_expr
    )

    # ============================================================
    # 46. In-plane stress output
    #
    # [sigma_xx, sigma_yy, sigma_xy]
    # ============================================================

    stress_2d_ufl = (
        tensor_to_voigt2(
            sigma_final
        )
    )

    stress_expr = fem.Expression(
        stress_2d_ufl,
        V_P1_voigt.element.interpolation_points
    )

    stress = fem.Function(
        V_P1_voigt
    )

    stress.interpolate(
        stress_expr
    )

    # ============================================================
    # 47. Total strain output
    #
    # [eps_xx, eps_yy, eps_xy]
    # ============================================================

    strain_2d_ufl = (
        tensor_to_voigt2(
            eps_total_final
        )
    )

    strain_expr = fem.Expression(
        strain_2d_ufl,
        V_P1_voigt.element.interpolation_points
    )

    strain_fn = fem.Function(
        V_P1_voigt
    )

    strain_fn.interpolate(
        strain_expr
    )

    # ============================================================
    # 48. Elastic strain output
    # ============================================================

    strain_el_2d_ufl = (
        tensor_to_voigt2(
            eps_el_final
        )
    )

    strain_el_expr = fem.Expression(
        strain_el_2d_ufl,
        V_P1_voigt.element.interpolation_points
    )

    strain_el_fn = fem.Function(
        V_P1_voigt
    )

    strain_el_fn.interpolate(
        strain_el_expr
    )

    # ============================================================
    # 49. Thermal strain output
    # ============================================================

    strain_th_2d_ufl = (
        tensor_to_voigt2(
            eps_th_final
        )
    )

    strain_th_expr = fem.Expression(
        strain_th_2d_ufl,
        V_P1_voigt.element.interpolation_points
    )

    strain_th_fn = fem.Function(
        V_P1_voigt
    )

    strain_th_fn.interpolate(
        strain_th_expr
    )

    # ============================================================
    # 50. Transformation strain output
    # ============================================================

    strain_tr_2d_ufl = (
        tensor_to_voigt2(
            eps_tr_final
        )
    )

    strain_tr_expr = fem.Expression(
        strain_tr_2d_ufl,
        V_P1_voigt.element.interpolation_points
    )

    strain_tr_fn = fem.Function(
        V_P1_voigt
    )

    strain_tr_fn.interpolate(
        strain_tr_expr
    )

    # ============================================================
    # 51. Plastic strain output
    #
    # Convert the quadrature history back to a nodal P1 field
    # ONLY for visualization/output.
    #
    # The actual constitutive history remains at quadrature
    # points.
    # ============================================================

    eps_p_2d_ufl = ufl.as_vector([
        eps_p_old[0],
        eps_p_old[1],
        eps_p_old[3]
    ])

    eps_p_expr = fem.Expression(
        eps_p_2d_ufl,
        V_P1_voigt.element.interpolation_points
    )

    eps_p_fn = fem.Function(
        V_P1_voigt
    )

    eps_p_fn.interpolate(
        eps_p_expr
    )

    # ============================================================
    # 52. Equivalent plastic strain output
    # ============================================================

    ep_eq_expr = fem.Expression(
        ep_eq_old,
        V_P1.element.interpolation_points
    )

    ep_eq_fn = fem.Function(
        V_P1
    )

    ep_eq_fn.interpolate(
        ep_eq_expr
    )

    # ============================================================
    # 53. Final yield function
    # ============================================================

    sigma_y_final = (
        sigma_y0
        + H * ep_eq_old
    )

    yield_final = (
        vm_final
        - sigma_y_final
    )

    yield_expr = fem.Expression(
        yield_final,
        V_P1.element.interpolation_points
    )

    yield_fn = fem.Function(
        V_P1
    )

    yield_fn.interpolate(
        yield_expr
    )

    # ============================================================
    # 54. Results
    # ============================================================

    res_dict = {

        "Displacement":
            u_final.x.array
            .reshape(-1, gdim)
            [fenics_to_xdmf_map],

        "Temperature":
            T_final.x.array
            [fenics_to_xdmf_map],

        "vonMises":
            vm_stress.x.array
            [fenics_to_xdmf_map],

        "Martensite":
            fM_final.x.array
            [fenics_to_xdmf_map],

        "Austenite":
            (
                1.0
                - fM_final.x.array
            )
            [fenics_to_xdmf_map],

        "Stress":
            stress.x.array.real
            .reshape(-1, 3)
            [fenics_to_xdmf_map],

        "Strain":
            strain_fn.x.array.real
            .reshape(-1, 3)
            [fenics_to_xdmf_map],

        "Strain_el":
            strain_el_fn.x.array.real
            .reshape(-1, 3)
            [fenics_to_xdmf_map],

        "Strain_th":
            strain_th_fn.x.array.real
            .reshape(-1, 3)
            [fenics_to_xdmf_map],

        "Strain_trip":
            strain_tr_fn.x.array.real
            .reshape(-1, 3)
            [fenics_to_xdmf_map],

        "PlasticStrain":
            eps_p_fn.x.array.real
            .reshape(-1, 3)
            [fenics_to_xdmf_map],

        "EquivalentPlasticStrain":
            ep_eq_fn.x.array.real
            [fenics_to_xdmf_map],

        "YieldFunction":
            yield_fn.x.array.real
            [fenics_to_xdmf_map],

        "Stress_hydrostatic":
            sh.x.array.real
            [fenics_to_xdmf_map],

        "Triaxiality":
            trax.x.array.real
            [fenics_to_xdmf_map],

        "Lode":
            lode.x.array.real
            [fenics_to_xdmf_map]
    }

    # ============================================================
    # 55. Write results
    # ============================================================

    adjustdatastream(
        res_dict,
        datapos="nodes",
        t_data=0.0
    )

    print(
        ">>> Quadrature-point J2 plasticity "
        "quench solve finished."
    )
