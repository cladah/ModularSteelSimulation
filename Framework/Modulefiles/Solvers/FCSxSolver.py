from mpi4py import MPI
import numpy as np
import sys
import os
from dolfinx import mesh, fem, io, nls
import ufl
import basix
from basix.ufl import element
from petsc4py import PETSc
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

# strain and stress
def eps(u):
    return ufl.sym(ufl.grad(u))
def voigt_strain(symm_tensor):
    return ufl.as_vector([symm_tensor[0, 0], symm_tensor[1, 1], 2.0 * symm_tensor[0, 1]])
def voigt_to_tensor(v):
    return ufl.as_tensor([[v[0], v[2] / 2.0], [v[2] / 2.0, v[1]]])

def sigma(u, minput):
    E = 200e9
    nu = 0.3
    mu = E / (2 * (1 + nu))
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    return lmbda * ufl.tr(eps(u)) * ufl.Identity(2) + 2.0 * mu * eps(u)
def vM(sig):
    return ufl.tr(sig)
def sig_h(sig):
    pass

def MartensiteForm(domain, minput, V_fM, funcs):
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
    Ms_data = readdatastream("KM_Ms_Martensite")
    beta_data = readdatastream("KM_b_Martensite")
    nodes = readdatastream("nodes")
    fMeq = 1.0

    # Inputfunctions
    uC = fem.Function(V_fM) # Amount of carbon in the material
    Ms = fem.Function(V_fM)
    beta = fem.Function(V_fM)
    Ms.name = "Martensite starttemp"
    beta.name = "KM beta variable"
    wtC_interp = LinearNDInterpolator(nodes, wtC_data, fill_value=0.0)
    Ms_interp = LinearNDInterpolator(nodes, Ms_data, fill_value=0.0)
    beta_interp = LinearNDInterpolator(nodes, beta_data, fill_value=0.0)

    # Interpolating to calculation points
    uC.interpolate(lambda x: wtC_interp(x[0], x[1]))
    Ms.interpolate(lambda x: Ms_interp(x[0], x[1]))
    beta.interpolate(lambda x: beta_interp(x[0], x[1]))

    fM_expr = fem.Expression(ufl.conditional(ufl.gt(Ms, T), fMeq - ufl.exp(-beta * (Ms - T)), 0.0),
                             V_fM.element.interpolation_points())

    #fM_cond = ufl.conditional(ufl.gt(Ms - T_out, 0.0), 1.0, 0.0)
    fM_cond = 1
    fM_prob = ((fM - fM_old) / dtime + ((fMeq - fM) * beta * (T - T_old) * (50*vM(sigma(u, minput))- 50*vM(sigma(u_old, minput))) / dtime)) * dfM
    fM_res = fM_cond * fM_prob * ufl.dx
    return fM_res, fM_expr

def DisplacementForm(domain, minput, V_u, funcs):
    wtC_data = readdatastream("Composition_C")
    fM = funcs["fM"]
    fM_old = funcs["fM_old"]
    dfM = funcs["dfM"]
    u = funcs["u"]
    u_old = funcs["u_old"]
    du = funcs["du"]

    # Material parameters
    E = 200e9
    nu = 0.3
    """
    # Deformation gradient
    F = ufl.variable(ufl.Identity(2) + ufl.grad(u))

    # Right Cauchy-Green tensor
    C = F.T * F

    # Invariants of deformation tensors
    I1 = ufl.tr(C)
    J = ufl.det(F)

    # Continuum equation
    psi = (mu / 2) * (I1 - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    """
    disp_res = ufl.inner(sigma(u, minput), eps(du)) * ufl.dx # Small strain formulation
    #disp_res = ufl.inner(ufl.grad(du, x, nu), P) * x[0] * ufl.dx
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
def FCSx4PB_Force(parent):
    print('Using FeniCSx solver for Mechanical FEM calculation')
    ginput = parent.ginput
    minput = parent.minput
    nodes = readdatastream("nodes")

    domain = readDatastreamMesh()
    domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0,0]), np.array([0.12, 0.012])], [60, 20], cell_type=mesh.CellType.quadrilateral)

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
    #bcs = []
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    T_res, s_old, s_expr = TemperatureForm(domain, minput, V_T, funcs)
    u_res = DisplacementForm(domain, minput, V_u, funcs)
    #Mart_exp = ufl.conditional(ufl.gt(Ms, T_out), 1 - ufl.exp(-beta * (Ms - T_out)), 0.0)

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
    #resultdict["Displacement"] = u_old.x.array
    adjustdatastream(resultdict, "nodes", current_t)


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

    grid.point_data["uy"] = uy_out.x.array.real  # attach solution

    #plotter = pv.Plotter()
    #plotter.view_xy()
    #plotter.add_mesh(grid, scalars="uy", cmap="viridis")
    #plotter.show(title="FEniCSx solution with PyVista")

    #adjustdatastream(data, "nodes", current_t)
def FCSx4PB_Quench(parent):
    print('Using FeniCSx solver for FEM quenching calculation')
    ginput = parent.ginput
    minput = parent.minput
    nodes = readdatastream("nodes")

    domain = readDatastreamMesh()
    domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0,0]), np.array([0.12, 0.012])], [30, 20], cell_type=mesh.CellType.triangle)

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

    U = fem.Function(V)
    dU = ufl.TrialFunction(V)

    # Defining reference values
    Tref = fem.Function(V_T)
    Tref.x.array[:] = 293.15

    # Defining current functions
    funcs = dict()

    (u, T, fM) = ufl.split(U)
    U_ = ufl.TestFunction(V)
    (du, dT, dfM) = ufl.split(U_)
    dU = ufl.TrialFunction(V)
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

    Tref = fem.Function(V_T)
    Tref.x.array[:] = 293.15

    # Defining bounderies
    fdim = domain.topology.dim - 1

    bottom_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: np.isclose(x[1], 0.0, atol=1e-8))
    top_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: np.isclose(x[1], 0.012, atol=1e-8))
    bottom_dofs = fem.locate_dofs_topological((V.sub(1), V_T), fdim, bottom_facets)
    top_dofs = fem.locate_dofs_topological((V.sub(1), V_T), fdim, top_facets) # Important with the tuple????
    bottom_Tbc = fem.dirichletbc(T_bc, bottom_dofs, V.sub(1))
    top_Tbc = fem.dirichletbc(T_bc, top_dofs, V.sub(1))

    bcs = [bottom_Tbc, top_Tbc]
    #bcs = []
    fM_res, fM_exp = MartensiteForm(domain, minput, V_fM, funcs)
    T_res, s_old, s_expr = TemperatureForm(domain, minput, V_T, funcs)

    #Mart_exp = ufl.conditional(ufl.gt(Ms, T_out), 1 - ufl.exp(-beta * (Ms - T_out)), 0.0)

    Res = T_res + fM_res
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

    current_t = 0.0

    #topology, cell_types, geometry = dolfinx.plot.vtk_mesh(V_T)
    #grid = pv.UnstructuredGrid(topology, cell_types, geometry)

    for loop_nr in range(50):
        print("Loop nr " + str(loop_nr + 1))
        current_t += dtime
        num_its, converged = solver.solve(U)
        assert converged
        U.x.scatter_forward()
        T_old.interpolate(fem.Expression(T, V_T.element.interpolation_points()))

        # Updating history variables
        s_old.interpolate(s_expr)
        # fM_old.interpolate(fM_expr)
        u_old.interpolate(fem.Expression(u, V_u.element.interpolation_points()))
        #fM_old.interpolate(fem.Expression(fM, V_fM.element.interpolation_points()))
        fM_old.interpolate(fM_exp)
        T_old.interpolate(fem.Expression(T, V_T.element.interpolation_points()))

        # Adding solutions to results
        u_out = U.sub(0).collapse()
        u_out.name = "Displacement"
        T_out = U.sub(1).collapse()
        T_out.name = "Temperature"
        fM_out = U.sub(2).collapse()
        fM_out.name = "Martensite"

        fM_out.interpolate(fM_exp)

    #grid.point_data["T"] = T_out.x.array.real  # attach solution
    #grid.point_data["fM"] = fM_out.x.array.real
    #plotter = pv.Plotter()
    #plotter.add_mesh(grid, scalars="fM", cmap="viridis")
    #plotter.show(title="FEniCSx solution with PyVista")

def Compression():
    pass
def Test(parent):
    import pyvista
    import pandas as pd
    import ufl
    import numpy as np
    import scipy
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from basix.ufl import element
    import basix

    from petsc4py import PETSc
    from mpi4py import MPI
    from dolfinx import fem, mesh, io, plot, default_scalar_type, log
    from dolfinx.fem.petsc import assemble_matrix, assemble_vector, create_vector, apply_lifting, set_bc
    from dolfinx import nls
    from dolfinx.nls.petsc import NewtonSolver
    log.set_log_level(log.LogLevel.ERROR)
    t0 = 0.0
    t = t0
    tend = 20.0
    t_steps = 50
    T0 = 1000.0 + 273.15
    Tsurf = 20.0 + 273.15

    k = 44.5  # W/(m*K)
    Cp = 550.0  # J/(kg*K)
    Cp_a = [[0, 800], [200, 600]]
    R = 8.314
    Cv = Cp - R
    alpha = 2.31e-5
    E = 210e9
    nu = 0.3
    sig0 = 1000.0e6  # yield strength in MPa
    Et = E / 100.0  # tangent modulus
    H = E * Et / (E - Et)  # hardening modulus

    # Deviatoric components
    lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
    mu = E / 2 / (1 + nu)
    kappa = alpha * (3 * lmbda + 2 * mu)
    kappa_fM = 0.004 * (3 * lmbda + 2 * mu)
    # Tref_value = 293.15 # Reference temperature
    Tref_value = T0
    fMeq = 0.9

    dim = 2

    th = 0.008  # m
    radius = 0.008  # m
    height = 0.001  # m
    rho = 7800  # kg/m3

    domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0,0]), np.array([0.007, 0.001])], [90, 30], cell_type=mesh.CellType.triangle)

    dt = fem.Constant(domain, tend / t_steps)

    # Defining elements spaces
    Vue = element("P", domain.basix_cell(), 2, shape=(2,))  # displacement finite element
    VTe = element("Lagrange", domain.basix_cell(), 1)  # temperature finite element
    VfMe = element("Lagrange", domain.basix_cell(), 1)  # martensite finite element
    V = fem.functionspace(domain, basix.ufl.mixed_element([Vue, VTe, VfMe]))

    V_u, _ = V.sub(0).collapse()
    V_ux, _ = V.sub(0).sub(0).collapse()  # used for Dirichlet BC
    V_uy, _ = V.sub(0).sub(1).collapse()  # used for Dirichlet BC
    V_T, _ = V.sub(1).collapse()  # used for Dirichlet BC
    V_fM, _ = V.sub(2).collapse()

    cell = domain.basix_cell()
    quadrature_rule = basix.make_quadrature(cell, 2)

    # basis_values = element.tabulate(0, quadrature_rule.points)
    # Storing reference information in fem.Functions
    Tref = fem.Function(V_T)
    Tref.x.array[:] = Tref_value

    # Martensite transformation data from csv files
    Ms_dat = pd.read_csv("Resultfiles/Martensite_Ms.csv", names=["Radius", "Ms"])
    Ms_dat["Ms"] = Ms_dat["Ms"]
    Ms_dat.loc[len(Ms_dat)] = {'radius': 0.0081, 'Ms': 300}
    beta_dat = pd.read_csv("Resultfiles/Martensite_beta.csv", names=["Radius", "Beta"])
    f_row = beta_dat.iloc[0].copy()
    f_row['Radius'] = f_row['Radius'] - beta_dat.iloc[1]["Radius"]
    e_row = beta_dat.iloc[-1].copy()
    e_row['Radius'] = e_row['Radius'] + beta_dat.iloc[-1]["Radius"] - beta_dat.iloc[-2]["Radius"]
    beta_dat = pd.concat([pd.DataFrame([f_row]), beta_dat, pd.DataFrame([e_row])], ignore_index=True)

    # Interpolation in 1D to radial coordinates
    Ms_func = scipy.interpolate.interp1d(Ms_dat["Radius"].round(8), Ms_dat["Ms"], fill_value=0.0)
    beta_func = scipy.interpolate.interp1d(beta_dat["Radius"].round(8), beta_dat["Beta"], fill_value=0.0)

    Ms = fem.Function(V_fM)
    Ms.interpolate(lambda x: Ms_func(np.sqrt(x[0] ** 2 + x[1] ** 2)))
    Ms.name = "Martensite starttemp"
    beta = fem.Function(V_fM)
    beta.name = "KM beta variable"
    beta.interpolate(lambda x: beta_func(np.sqrt(x[0] ** 2 + x[1] ** 2)))

    # Defining result functions
    x = ufl.SpatialCoordinate(domain)
    T_out = fem.Function(V_T)
    fM_out = fem.Function(V_fM)
    ux_out = fem.Function(V)
    uy_out = fem.Function(V)

    Mart_exp = ufl.conditional(ufl.gt(Ms, T_out), 1 - ufl.exp(-beta * (Ms - T_out)), 0.0)
    # Locating Boundaries
    fdim = domain.topology.dim - 1

    surface_facets = mesh.locate_entities_boundary(domain, fdim,
                                                   lambda x: np.isclose(np.sqrt(x[0] ** 2 + x[1] ** 2), radius))
    bottom_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: np.isclose(x[1], 0.0))
    left_facets = mesh.locate_entities_boundary(domain, fdim, lambda x: np.isclose(x[0], 0.0))

    surface_T = fem.locate_dofs_topological((V.sub(1), V_T), fdim, surface_facets)
    bottom_uy = fem.locate_dofs_topological((V.sub(0).sub(1), V_uy), fdim, bottom_facets)
    left_ux = fem.locate_dofs_topological((V.sub(0).sub(0), V_ux), fdim, left_facets)

