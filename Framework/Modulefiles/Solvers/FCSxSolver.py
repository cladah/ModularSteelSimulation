from mpi4py import MPI
import numpy as np
import sys
import os
from dolfinx import mesh, fem, io, nls
import ufl
import basix
import adios4dolfinx
import dolfinx
from dolfinx.fem import Function, FunctionSpace, Constant, dirichletbc, form, locate_dofs_topological, assemble_scalar, assemble_vector
from Framework.HelpFile import readresultfile, read_modinput, read_geninput
from Framework.Datastream_file import readdatastream, readDatastreamMesh, adjustdatastream
#from ResultReading import getnames_results
import meshio
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from scipy.interpolate import LinearNDInterpolator
import pyvista as pv
from dolfinx.io import XDMFFile
from dolfinx import fem, mesh, plot, default_scalar_type
from basix.ufl import element, mixed_element
from petsc4py import PETSc

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

def load_field_to_function(domain, filename, field_name):
    """
    Loads a specific data field from the XDMF into a FEniCSx Function.
    """
    # Define a scalar FunctionSpace (Lagrange Degree 1 for node data)
    V = dolfinx.fem.functionspace(domain, ("Lagrange", 1))
    u = dolfinx.fem.Function(V)
    u.name = field_name

    with XDMFFile(domain.comm, filename, "r") as xdmf:
        try:
            # Modern versions use read_checkpoint for Function data
            xdmf.read_checkpoint(u, field_name, 0)
        except AttributeError:
            # Fallback for slightly different dev versions
            # Some versions use a mesh-specific read
            xdmf.read_mesh(name=field_name)

    return u

# strain and stress

def voigt_strain(symm_tensor):
    return ufl.as_vector([symm_tensor[0, 0], symm_tensor[1, 1], 2.0 * symm_tensor[0, 1]])
def voigt_to_tensor(v):
    return ufl.as_tensor([[v[0], v[2] / 2.0], [v[2] / 2.0, v[1]]])


def eps(u):
    return ufl.sym(ufl.grad(u))

def sigma(u, minput):
    E, nu = 210e9, 0.3
    mu = E / (2 * (1 + nu))
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lmbda * ufl.div(u) * ufl.Identity(len(u)) + 2 * mu * eps(u)

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
    #Ms = load_field_to_function(domain, "Datastream.xdmf", "KM_Ms_Martensite")
    #beta = load_field_to_function(domain, "Datastream.xdmf", "KM_b_Martensite")
    Ms = 350
    beta = 0.01
    fMeq = 1.0

    # Koistinen-Marburger Expression for current state
    fM_target = ufl.conditional(ufl.gt(Ms, T),
                                fMeq * (1.0 - ufl.exp(-beta * (Ms - T + 1e-7 * triax(sigma(u, minput))*vM(sigma(u, minput))))),
                                0.0)

    fM_expr = fem.Expression(fM_target, V_fM.element.interpolation_points())

    # Residual: Standard L2 Projection of the algebraic KM equation
    fM_res = (fM - fM_target) * dfM * ufl.dx

    return fM_res, fM_expr

def MartensiteForm_old(domain, minput, V_fM, funcs):
    fM = funcs["fM"]
    fM_old = funcs["fM_old"]
    dfM = funcs["dfM"]
    T = funcs["T"]
    T_old = funcs["T_old"]
    dT = funcs["dT"]
    u = funcs["u"]
    u_old = funcs["u_old"]
    du = funcs["du"]
    dtime = 0.01
    # Reading datastream data
    wtC_data = readdatastream("Composition_C")
    #Ms_data = readdatastream("KM_Ms_Martensite")
    #beta_data = readdatastream("KM_b_Martensite")
    nodes = readdatastream("nodes")
    fMeq = 1.0

    # Inputfunctions
    uC = fem.Function(V_fM) # Amount of carbon in the material
    Ms = fem.Function(V_fM)
    beta = fem.Function(V_fM)
    Ms.name = "Martensite starttemp"
    beta.name = "KM beta variable"
    wtC_interp = LinearNDInterpolator(nodes, wtC_data, fill_value=0.0)
    #Ms_interp = LinearNDInterpolator(nodes, Ms_data, fill_value=0.0)
    #beta_interp = LinearNDInterpolator(nodes, beta_data, fill_value=0.0)

    # Interpolating to calculation points
    uC.interpolate(lambda x: wtC_interp(x[0], x[1]))
    #Ms.interpolate(lambda x: Ms_interp(x[0], x[1]))
    #beta.interpolate(lambda x: beta_interp(x[0], x[1]))

    fM_expr = fem.Expression(ufl.conditional(ufl.gt(Ms, T), fMeq - ufl.exp(-beta * (Ms - T)), 0.0),
                             V_fM.element.interpolation_points())
    fM_expr = fem.Expression(fMeq - ufl.exp(-beta * (0 - T)), V_fM.element.interpolation_points())

    #fM_cond = ufl.conditional(ufl.gt(Ms - T_out, 0.0), 1.0, 0.0)
    fM_cond = 1
    #fM_prob = ((fM - fM_old) / dtime + ((fMeq - fM) * beta * (T - T_old) * (50*vM(sigma(u, minput))- 50*vM(sigma(u_old, minput))) / dtime)) * dfM
    fM_prob = (fMeq - fM) * (50*vM(sigma(u, minput)))/ dtime*dfM
    fM_res = fM_cond * fM_prob * ufl.dx
    return fM_res, fM_expr

def DisplacementForm(domain, minput, V_u, funcs):
    wtC_data = readdatastream("Composition_C")
    fM = funcs["fM"]
    fM_old = funcs["fM_old"]
    dfM = funcs["dfM"]
    u = funcs["u"]
    du = funcs["du"]

    disp_res = ufl.inner(sigma(u, minput), eps(du)) * ufl.dx
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

def TemperatureForm(domain, minput, V_T, funcs):
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

    # Entropyy
    s = Cv / Tref * T
    s_expr = fem.Expression(s, V_T.element.interpolation_points())
    s_old = fem.Function(V_T)
    s_old.x.array[:] = (T0 * Cv) / Tref
    # Write a script that gets top and bottom values in y
    surface_facets = mesh.locate_entities_boundary(domain, 1,
                                                   lambda x: np.logical_or(np.isclose(x[1], 0.0), np.isclose(x[1], 0.012)))
    facet_tag = mesh.meshtags(domain, 1, surface_facets, np.full(len(surface_facets), 1, dtype=np.int32))
    surf_ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag)
    therm_res = (
                        rho * Tref * (s - s_old) / dtime * dT - ufl.dot(-k * ufl.grad(T), ufl.grad(dT))
                ) * ufl.dx - 20000.0 * ufl.inner((Tsurf-T), dT) * surf_ds(1)

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
    print('>>> Using FEniCSx solver: 4-Point Bend (Displacement/Force Control)')
    minput = parent.minput

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
    T_old.x.array[:] = 293.15
    fM_old = fem.Function(V_fM)

    funcs = {"u": u, "du": du, "u_old": u_old, "T": T, "dT": dT,
             "T_old": T_old, "fM": fM, "dfM": dfM, "fM_old": fM_old, "dtime": dtime}

    # 5. Boundary Conditions (4-Point Bending)
    def b1_BC(x): return np.logical_and(np.isclose(x[0], 0.0), np.isclose(x[1], 0.0))

    def t1_BC(x): return np.logical_and(np.isclose(x[0], 0.04), np.isclose(x[1], 0.012))

    def sym_x_BC(x): return np.isclose(x[0], 0.06)

    def sym_z_BC(x): return np.isclose(x[2], 0.00)

    # Dirichlet setup
    uT_val = fem.Function(V_uy)
    uT_val.x.array[:] = 0.0
    ux_sym_val = fem.Function(V_ux)
    ux_sym_val.x.array[:] = 0.0
    uz_sym_val = fem.Function(V_uz)
    uz_sym_val.x.array[:] = 0.0
    u_zero_func = fem.Function(V_u)
    u_zero_func.x.array[:] = 0.0

    bcs = [
        fem.dirichletbc(ux_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), b1_BC), V.sub(0).sub(1)),
        fem.dirichletbc(uT_val, fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), t1_BC), V.sub(0).sub(1)),
        fem.dirichletbc(ux_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(0), V_ux), sym_x_BC), V.sub(0).sub(0)),
        fem.dirichletbc(uz_sym_val, fem.locate_dofs_geometrical((V.sub(0).sub(2), V_uz), sym_z_BC), V.sub(0).sub(2)),
    ]

    # 6. Variational Forms
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    #T_res, s_old, s_expr = TemperatureForm(domain, minput, V_T, funcs)
    T_dummy = (T - T_old) * dT * ufl.dx
    fM_dummy = (fM - fM_old) * dfM * ufl.dx
    u_res = DisplacementForm(domain, minput, V_u, funcs)

    Res = u_res + T_dummy + fM_res
    Jac = ufl.derivative(Res, U, dU)

    # 7. Solver Setup
    problem = NonlinearProblem(Res, U, bcs=bcs, J=Jac)
    solver = NewtonSolver(domain.comm, problem)
    solver.convergence_criterion = "incremental"
    solver.max_it = 20
    solver.report = True

    # Configure PETSc for direct solve (more stable for these point BCs)
    ksp = solver.krylov_solver
    opts = PETSc.Options()
    prefix = ksp.getOptionsPrefix()
    opts[f"{prefix}ksp_monitor"] = None
    opts[f"{prefix}ksp_type"] = "preonly"
    opts[f"{prefix}pc_type"] = "lu"
    opts[f"{prefix}pc_factor_mat_solver_type"] = "mumps"
    ksp.setFromOptions()

    # 8. Time Loop
    print(">>> Starting Solve...")
    for step in range(minput["loadsteps"]):
        print(f"Step nr {step + 1}")
        uT_val.x.array[:] = -0.0005*(step+1)

        num_its, converged = solver.solve(U)
        U.x.scatter_forward()

        # Update History
        T_old.interpolate(U.sub(1))
        u_old.interpolate(U.sub(0))
        #s_old.interpolate(s_expr)
        #fM_old.interpolate(fem.Expression(fM_exp, V_fM.element.interpolation_points()))
    # 9. Extract Results for Datastream
    # Collapse to extract arrays
    u_final = U.sub(0).collapse()
    T_final = U.sub(1).collapse()
    fM_final = U.sub(2).collapse()

    vm_expr_ufl = sigma_von_mises(U.sub(0), minput)
    vm_expr = fem.Expression(vm_expr_ufl, V_P1.element.interpolation_points())
    vm_stress = fem.Function(V_P1)
    vm_stress.interpolate(vm_expr)

    lode_expr_ufl = lode_func(sigma(U.sub(0),minput))
    lode_expr = fem.Expression(lode_expr_ufl, V_P1.element.interpolation_points())
    lode = fem.Function(V_P1)
    lode.interpolate(lode_expr)

    trax_expr_ufl = triax(sigma(U.sub(0), minput))
    trax_expr = fem.Expression(trax_expr_ufl, V_P1.element.interpolation_points())
    trax = fem.Function(V_P1)
    trax.interpolate(trax_expr)

    sh_expr_ufl = sig_h(sigma(U.sub(0), minput))
    sh_expr = fem.Expression(sh_expr_ufl, V_P1.element.interpolation_points())
    sh = fem.Function(V_P1)
    sh.interpolate(sh_expr)

    s_expr_ufl = sigma(U.sub(0), minput)
    sv_expr_ufl = tensor_to_voigt(s_expr_ufl, dim=gdim)
    stress_expr = fem.Expression(sv_expr_ufl, V_P1_voigt.element.interpolation_points())
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
    plot_4PB_results(U, minput)
def FCSx4PB_Force_old(parent):
    print('Using FeniCSx solver for Mechanical FEM calculation. Force controlled.')
    ginput = parent.ginput
    minput = parent.minput
    nodes = readdatastream("nodes")

    domain = readDatastreamMesh()
    domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0,0]), np.array([0.12, 0.012])], [60, 20], cell_type=mesh.CellType.quadrilateral)
    with io.XDMFFile(MPI.COMM_WORLD, "Resultfiles/Mesh.xdmf", "r") as file:
        domain = file.read_mesh(name="Grid")
    print(np.max(domain.geometry.x[:,1]))
    xdmfpos = readdatastream('nodes')

    matrix1_expanded = xdmfpos[:, np.newaxis, :]
    differences = matrix1_expanded - domain.geometry.x[:, :2]
    distances = np.linalg.norm(differences, axis=2)
    best_match_indices = np.argmin(distances, axis=1)

    dtime = fem.Constant(domain, minput["quenchtime"]/minput["quench_steps"])
    tend = fem.Constant(domain, float(minput["quenchtime"]))

    Vue = element("P", domain.basix_cell(), 2, shape=(2,))  # displacement finite element
    VTe = element("Lagrange", domain.basix_cell(), 1)  # temperature finite element
    VfMe = element("Lagrange", domain.basix_cell(), 1)  # martensite finite element
    V = fem.functionspace(domain, basix.ufl.mixed_element([Vue, VTe, VfMe]))
    V_u, _ = V.sub(0).collapse()
    V_ux, _ = V.sub(0).sub(0).collapse()  # used for Dirichlet BC
    V_uy, _ = V.sub(0).sub(1).collapse()  # used for Dirichlet BC
    V_T, _ = V.sub(1).collapse()  # used for Dirichlet BC
    V_fM, _ = V.sub(2).collapse()

    V_post = fem.functionspace(domain, element("P", domain.basix_cell(), 1, shape=(2,)))

    # Defining reference values
    Tref = fem.Function(V_T)
    Tref.x.array[:] = 293.15

    # Defining current functions
    funcs = dict()

    U = fem.Function(V)
    dU = ufl.TrialFunction(V)
    (u, T, fM) = ufl.split(U)
    U_ = ufl.TestFunction(V)
    (du, dT, dfM) = ufl.split(U_)

    # Displacement function
    funcs["u"] = u
    funcs["du"] = du
    # Temperature
    funcs["T"] = T
    funcs["dT"] = dT
    # Martensite function
    funcs["fM"] = fM
    funcs["dfM"] = dfM
    # Adding timestep to function dict
    funcs["dtime"] = dtime


    # Defining history functions
    fM_old = fem.Function(V_fM)
    fM_old.x.array[:] = 0.0
    funcs["fM_old"] = fM_old
    T_old = fem.Function(V_T)
    T_old.x.array[:] = 1113.15
    funcs["T_old"] = T_old
    u_old = fem.Function(V_u)
    u_old.x.array[:] = 0.0
    funcs["u_old"] = u_old

    # Defining result functions
    u_out = fem.Function(V_u)
    u_out.name = "Displacement"
    T_out = fem.Function(V_T)
    T_out.name = "Temperature"
    fM_out = fem.Function(V_fM)
    fM_out.name = "Martensite"

    # Defining boundary values
    T_bc = fem.Function(V_T)
    T_bc.x.array[:] = 293.15
    u_bc = fem.Function(V_u)
    u_bc.x.array[:] = 0.0
    uT_bc = fem.Function(V_uy)
    uT_bc.x.array[:] = 0.0

    Tref = fem.Function(V_T)
    Tref.x.array[:] = 293.15

    # Defining bounderies
    fdim = domain.topology.dim - 1


    def b1_BC(x):
        return np.logical_and(np.isclose(x[0], 0.0, atol=1e-8),np.isclose(x[1], 0.0, atol=1e-8))
    b1_dofs = fem.locate_dofs_geometrical((V.sub(0), V_u), b1_BC)
    b1_bc = fem.dirichletbc(u_bc, b1_dofs, V.sub(0))

    def b2_BC(x):
        return np.logical_and(np.isclose(x[0], 0.12, atol=1e-8), np.isclose(x[1], 0.0, atol=1e-8))
    b2_dofs = fem.locate_dofs_geometrical((V.sub(0), V_u), b2_BC)
    b2_bc = fem.dirichletbc(u_bc, b2_dofs, V.sub(0))

    def t1_BC(x):
        return np.logical_and(np.isclose(x[0], 0.04, atol=1e-8),np.isclose(x[1], 0.012, atol=1e-8))
    t1_dofs = fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), t1_BC)
    t1_bc = fem.dirichletbc(uT_bc, t1_dofs, V.sub(0).sub(1))

    def t2_BC(x):
        return np.logical_and(np.isclose(x[0], 0.08, atol=1e-8), np.isclose(x[1], 0.012, atol=1e-8))
    t2_dofs = fem.locate_dofs_geometrical((V.sub(0).sub(1), V_uy), t2_BC)
    t2_bc = fem.dirichletbc(uT_bc, t2_dofs, V.sub(0).sub(1))

    bcs = [b1_bc, b2_bc, t1_bc, t2_bc]
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    T_res, s_old, s_expr = TemperatureForm(domain, minput, V_T, funcs)
    u_res = DisplacementForm(domain, minput, V_u, funcs)
    #Mart_exp = ufl.conditional(ufl.gt(Ms, T_out), 1 - ufl.exp(-beta * (Ms - T_out)), 0.0)

    Res = u_res+T_res+fM_res
    Res = u_res
    Jac = ufl.derivative(Res, U, dU)

    problem = fem.petsc.NonlinearProblem(Res, U, bcs=bcs, J=Jac)

    print("Starting the solve")
    solver = NewtonSolver(domain.comm, problem)
    solver.convergence_criterion = "incremental"
    solver.report = True
    solver.max_it = 30
    ksp = solver.krylov_solver

    opts = PETSc.Options()

    option_prefix = ksp.getOptionsPrefix()
    opts[f"{option_prefix}ksp_type"] = "gmres"
    opts[f"{option_prefix}ksp_rtol"] = 1e-12
    opts[f"{option_prefix}ksp_atol"] = 1e-12
    opts[f"{option_prefix}pc_type"] = "hypre"
    opts[f"{option_prefix}pc_hypre_type"] = "boomeramg"
    opts[f"{option_prefix}pc_hypre_boomeramg_max_iter"] = 5
    opts[f"{option_prefix}pc_hypre_boomeramg_cycle_type"] = "v"
    ksp.setFromOptions()

    U.sub(0).interpolate(u_old)
    U.sub(1).interpolate(T_old)
    U.sub(2).interpolate(fM_old)
    U.x.scatter_forward()

    num_its, converged = solver.solve(U)
    U.sub(2).interpolate(fM_old)
    assert converged
    U.x.scatter_forward()
    T_old.interpolate(fem.Expression(T, V_T.element.interpolation_points()))

    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(V_T)
    grid = pv.UnstructuredGrid(topology, cell_types, geometry)
    grid.point_data["T"] = T_old.x.array.real  # attach solution

    U.sub(0).interpolate(u_old)
    U.sub(1).interpolate(T_old)
    U.sub(2).interpolate(fM_old)
    U.x.scatter_forward()

    current_t = 0.0

    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(V_uy)
    grid = pv.UnstructuredGrid(topology, cell_types, geometry)
    resultdict = dict()
    resultdict["T"] = T_old.x.array
    resultdict["Martensite"] = fM_old.x.array
    resultdict["Austenite"] = np.ones(len(fM_old.x.array)) - fM_old.x.array
    resultdict["Ferrite"] = np.zeros(len(fM_old.x.array))
    resultdict["Bainite"] = np.zeros(len(fM_old.x.array))
    resultdict["Pearlite"] = np.zeros(len(fM_old.x.array))
    u_nodes = fem.Function(V_post)


    for loop_nr in range(2):
        print("Loop nr " + str(loop_nr + 1))
        current_t += dtime
        num_its, converged = solver.solve(U)
        assert converged
        U.x.scatter_forward()
        T_old.interpolate(fem.Expression(T, V_T.element.interpolation_points()))
        uT_bc.x.array[:] = uT_bc.x.array[:] + 0.0001
        # Updating history variables
        s_old.interpolate(s_expr)
        # fM_old.interpolate(fM_expr)
        u_old.interpolate(fem.Expression(u, V_u.element.interpolation_points()))
        # fM_old.interpolate(fem.Expression(fM, V_fM.element.interpolation_points()))
        fM_old.interpolate(fM_exp)
        T_old.interpolate(fem.Expression(T, V_T.element.interpolation_points()))

        # Adding solutions to results
        u_out = U.sub(0).collapse()
        u_out.name = "Displacement"
        uy_out = U.sub(0).sub(1).collapse()
        uy_out.name = "Displacement in y"
        T_out = U.sub(1).collapse()
        T_out.name = "Temperature"
        fM_out = U.sub(2).collapse()
        fM_out.name = "Martensite"
        fM_out.interpolate(fM_exp)
        u_nodes.interpolate(u_out)
    grid.point_data["uy"] = uy_out.x.array.real  # attach solution

    resultdict["Displacement"] = u_nodes.x.array[best_match_indices]
    adjustdatastream(resultdict, datapos="nodes", t_data=0.0)
    print("FeniCSx testing")
    print(np.max(u_nodes.x.array))
    print(np.shape(fM_out.x.array))
    print(np.shape(u_nodes.x.array))
    print(np.shape(u_out.x.array))

    plotter = pv.Plotter()
    plotter.view_xy()
    plotter.add_mesh(grid, scalars="uy", cmap="viridis")
    plotter.show(title="FEniCSx solution with PyVista")

    #adjustdatastream(data, "nodes", current_t)

def FCSx4PB_Quench(parent):
    print('>>> Using FEniCSx solver for FEM quenching calculation')
    ginput = parent.ginput
    minput = parent.minput

    # 1. Load Domain
    domain = load_datastream_to_fenicsx("Datastream.xdmf")
    print("Domain dimensions!")
    print(domain.geometry.dim)
    print("       ")
    # 2. Setup Time Stepping
    dt_val = minput["quenchtime"] / minput["quench_steps"]
    dtime = fem.Constant(domain, default_scalar_type(dt_val))

    # 3. Define Mixed Function Space
    # Displacement (Vector), Temperature (Scalar), Martensite (Scalar)
    P2 = element("Lagrange", domain.basix_cell(), 2, shape=(domain.geometry.dim,))
    P1 = element("Lagrange", domain.basix_cell(), 1)

    V_el = mixed_element([P2, P1, P1])
    V = fem.functionspace(domain, V_el)

    # Sub-spaces for collapsing/interpolation
    V_u, _ = V.sub(0).collapse()
    V_T, _ = V.sub(1).collapse()
    V_fM, _ = V.sub(2).collapse()

    # 4. Functions
    U = fem.Function(V)
    u, T, fM = ufl.split(U)
    U_test = ufl.TestFunction(V)
    du, dT, dfM = ufl.split(U_test)
    dU_trial = ufl.TrialFunction(V)

    # History Functions (Stored at previous timestep)
    T_old = fem.Function(V_T)
    T_old.x.array[:] = 1113.15  # initial temperature

    fM_old = fem.Function(V_fM)
    fM_old.x.array[:] = 0.0

    u_old = fem.Function(V_u)
    u_old.x.array[:] = 0.0

    # Pack into dict for your Form functions
    funcs = {
        "u": u, "du": du, "u_old": u_old,
        "T": T, "dT": dT, "T_old": T_old,
        "fM": fM, "dfM": dfM, "fM_old": fM_old,
        "dtime": dtime
    }

    # 5. Boundary Conditions
    fdim = domain.topology.dim - 1
    # Assuming height of 0.012 from your code
    top_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: np.isclose(x[1], 0.012))
    bottom_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: np.isclose(x[1], 0.0))

    T_bc_val = fem.Function(V_T)
    T_bc_val.x.array[:] = 293.15

    dofs_top = fem.locate_dofs_topological((V.sub(1), V_T), fdim, top_facets)
    dofs_bottom = fem.locate_dofs_topological((V.sub(1), V_T), fdim, bottom_facets)

    bcs = [
        fem.dirichletbc(T_bc_val, dofs_top, V.sub(1)),
        fem.dirichletbc(T_bc_val, dofs_bottom, V.sub(1))
    ]

    # 6. Forms (Assuming these exist in your Framework)
    #fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    #T_res, s_old, s_expr = TemperatureForm(domain, minput, V_T, funcs)
    u_res = DisplacementForm(domain, minput, V_u, funcs)
    Res = u_res
    Jac = ufl.derivative(Res, U, dU_trial)

    # 7. Solver Setup
    problem = fem.petsc.NonlinearProblem(Res, U, bcs=bcs, J=Jac)
    solver = NewtonSolver(domain.comm, problem)
    solver.convergence_criterion = "incremental"
    solver.max_it = 50

    # KSP/PETSc Options
    ksp = solver.krylov_solver
    opts = PETSc.Options()
    prefix = ksp.getOptionsPrefix()
    opts[f"{prefix}ksp_type"] = "preonly"
    opts[f"{prefix}pc_type"] = "lu"  # Direct solver usually more stable for mixed forms
    opts[f"{prefix}pc_factor_mat_solver_type"] = "mumps"
    ksp.setFromOptions()

    # Initial Condition Interpolation
    U.sub(1).interpolate(T_old)
    U.x.scatter_forward()

    # 8. Time Loop
    print(">>> Starting Solve...")
    for step in range(int(minput["quench_steps"])):
        print(f"Time Step {step + 1}/{minput['quench_steps']}")

        num_its, converged = solver.solve(U)
        U.x.scatter_forward()

        # Update History
        T_old.interpolate(U.sub(1))
        fM_old.interpolate(fem.Expression(fM_exp, V_fM.element.interpolation_points()))
        u_old.interpolate(U.sub(0))
        s_old.interpolate(s_expr)

    # 9. Plotting Result
    plot_quenching_results(U)

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
    print(vm_stress.x.array.real)
    print(np.max(vm_stress.x.array.real))
    print(np.min(vm_stress.x.array.real))
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