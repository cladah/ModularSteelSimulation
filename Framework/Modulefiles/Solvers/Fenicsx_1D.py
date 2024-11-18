
import numpy as np
import ufl
import gmsh
from mpi4py import MPI
from dolfinx import fem, io, nls, default_scalar_type, log
import dolfinx.fem.petsc
from dolfinx.fem.petsc import LinearProblem
import dolfinx.nls.petsc
from dolfinx.mesh import create_box, CellType
from ufl import (
    as_matrix,
    dot,
    cos,
    sin,
    SpatialCoordinate,
    Identity,
    grad,
    ln,
    tr,
    det,
    variable,
    derivative,
    TestFunction,
    TrialFunction,
)
from basix.ufl import element
from dolfinx.io import XDMFFile, gmshio
import gmsh  # type: ignore
import meshio

def createxdmf():
    filename = "FCSX.xdmf"
    mesh = meshio.read("FeniCSx/FCSX_Mesh.msh")
    print(mesh.get_cells_type)
    with meshio.xdmf.TimeSeriesWriter(filename) as writer:
        writer.write_points_cells(points=mesh.points, cells={"line": mesh.get_cells_type("line")})
        writer.write_data(0, cell_data={})

def adjustxdmf(data, datapos="nodes", t_data=0.0):
    # Adding data to xdmf file
    try:
        if type(data) == dict:
            pass
        else:
            raise TypeError
    except TypeError:
        raise KeyError("datatype for datastream storage should be a dict")

    try:
        with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
            points, cells = reader.read_points_cells()
            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
    except meshio._exceptions.ReadError:
        raise KeyError("No meshgrid for timeseries")
    if t_data not in t_list:
        t_list.append(t_data)
        pd_list.append(data)
        cd_list.append({})
    else:
        indx = t_list.index(t_data)
        for key in data.keys():
            pd_list[indx][key] = data[key]
    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(points, cells)
        for i in range(len(t_list)):
            writer.write_data(t=t_list[i], point_data=pd_list[i], cell_data=cd_list[i])
    #print("Added " + str(dataname) + " to datastream")
    return

def gmsh1D():

    print('Remeshing 1D')

    data = dict()
    data["Geometry"] = {"nodes": 5, "meshscaling": 1, "radius": 0.008}

    r = data['Geometry']['radius']
    tmpgeo = [data['Geometry']['meshscaling'] ** i for i in range(data['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(data['Geometry']['nodes'] - 1)])
    rnodes = r * tmpgeo / np.max(tmpgeo)
    lc = rnodes[-1] - rnodes[-2]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()

    gmsh.clear()
    gmsh.model.add("Line")
    gdim = 3
    gmsh.option.setNumber("Geometry.Tolerance", 1.E-6)

    # Adding three points
    gmsh.model.occ.addPoint(0, 0, 0, 1)
    gmsh.model.occ.addPoint(r, 0, 0, lc, 2)
    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(1, [1], 1, 'Radius')

    # Geometric scaling of mesh
    gmsh.model.mesh.set_transfinite_curve(1, data['Geometry']['nodes'], 'Progression',
                                          data['Geometry']['meshscaling'])

    # Generating mesh
    gmsh.model.mesh.generate(gdim)
    #gmsh.model.mesh.setOrder(1)
    # ----------------------
    print(*gmsh.logger.get(), sep="\n")
    gmsh.write("FeniCSx/FCSX_Mesh.msh")
    print("Created 1D mesh for FenicsX")

    return gmsh.model

def gmsh1D_Really1D():

    print('Remeshing 1D')

    data = dict()
    data["Geometry"] = {"nodes": 5, "meshscaling": 1, "radius": 0.008}

    r = data['Geometry']['radius']
    tmpgeo = [data['Geometry']['meshscaling'] ** i for i in range(data['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(data['Geometry']['nodes'] - 1)])
    rnodes = r * tmpgeo / np.max(tmpgeo)
    lc = rnodes[-1] - rnodes[-2]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()

    gmsh.clear()
    gmsh.model.add("Line")
    gdim = 1
    gmsh.option.setNumber("Geometry.Tolerance", 1.E-6)

    # Adding three points
    gmsh.model.occ.addPoint(0, 0, 0, 1)
    gmsh.model.occ.addPoint(r, 0, 0, lc, 2)
    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(1, [1], 1, 'Radius')

    # Geometric scaling of mesh
    gmsh.model.mesh.set_transfinite_curve(1, data['Geometry']['nodes'], 'Progression',
                                          data['Geometry']['meshscaling'])

    # Generating mesh
    gmsh.model.mesh.generate(gdim)
    #gmsh.model.mesh.setOrder(1)
    # ----------------------
    print(*gmsh.logger.get(), sep="\n")
    gmsh.write("FeniCSx/FCSX_Mesh.msh")
    print("Created 1D mesh for FenicsX")

    return gmsh.model

def FNXTest():
    # +
    from dolfinx import log, default_scalar_type
    from dolfinx.fem.petsc import NonlinearProblem
    from dolfinx.nls.petsc import NewtonSolver
    import numpy as np
    import ufl

    from mpi4py import MPI
    from dolfinx import fem, mesh, plot
    import basix

    # Creating calculation grid
    gmshmodel = gmsh1D()
    createxdmf()

    gdim = 3
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)
    V = fem.functionspace(domain, ("Lagrange", 2, (gdim,)))
    VT = fem.functionspace(domain, ("Lagrange", 2, (gdim, gdim)))
    Vscal = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    x = SpatialCoordinate(domain)

    # Elasticity parameters
    E = default_scalar_type(200.0e9)
    nu = default_scalar_type(0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
    alpha = default_scalar_type(2.2e-05)
    T0 = default_scalar_type(1073.15)
    TRef = default_scalar_type(1073.15)


    # -

    # We then create a function for supplying the boundary condition on the left side, which is fixed.

    u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)

    # To apply the boundary condition, we identity the dofs located on the facets marked by the `MeshTag`.
    def center(x):
        return np.isclose(x[0], 0)

    def circ(x):
        return np.isclose(x[0], 0.008)

    center_dofs = fem.locate_dofs_geometrical(V, center)
    circ_dofs = fem.locate_dofs_geometrical(V, circ)
    center_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
    # V.sub(0)
    #bcs = [fem.dirichletbc(u_bc, center_dofs, V), fem.dirichletbc(u_bc, circ_dofs, V)]
    bcs = [fem.dirichletbc(u_bc, center_dofs, V)]

    # Next, we define the body force on the reference configuration (`B`), and nominal (first Piola-Kirchhoff) traction (`T`).
    BodF = fem.Constant(domain, default_scalar_type((0, 0, 0)))
    Trac = fem.Constant(domain, default_scalar_type((0, 0, 0)))

    # Define the test and solution functions on the space $V$

    u = fem.Function(V)
    u.name = "Displacement"
    v = ufl.TestFunction(V)

    Temp = fem.Function(Vscal)
    Temp.name = "Temperature"
    dTemp = ufl.TestFunction(Vscal)

    fM = fem.Function(Vscal)
    fM.name = "Martensite"
    dfM = ufl.TestFunction(Vscal)

    fB = fem.Function(Vscal)
    fB.name = "Bainate"
    dfB = ufl.TestFunction(Vscal)

    fP = fem.Function(Vscal)
    fP.name = "Pearlite"
    dfP = ufl.TestFunction(Vscal)

    fA = fem.Function(Vscal)
    fA.name = "Austenite"
    dfA = ufl.TestFunction(Vscal)




    # Define kinematic quantities used in the problem

    # Spatial dimension
    d = len(u)

    # Identity tensor
    I = ufl.variable(ufl.Identity(d))

    # Deformation gradient
    F = ufl.variable(I + ufl.grad(u))
    gradxyz = ufl.grad(u)
    F_cyl = ufl.variable(I + ufl.as_tensor([[gradxyz[0, 0], (gradxyz[0, 1] - u[1])/x[0], gradxyz[0, 2]],
                                            [gradxyz[1, 0], (gradxyz[1, 1] + u[0])/x[0], gradxyz[1, 2]/x[0]],
                                            [gradxyz[2, 0], gradxyz[2, 1]/x[0], gradxyz[2, 2]]]))
    F_cyl = ufl.variable(I + ufl.as_tensor([[u[0]*x[0], 0, 0],
                                            [0, u[0] / x[0], 0],
                                            [0, 0, 0]]))
    #F_cyl = ufl.variable(I + ufl.as_tensor([[gradxyz[0, 0], (gradxyz[0, 1] - u[1]) / x[0], 0],
    #                                        [gradxyz[1, 0], (gradxyz[1, 1] + u[0]) / x[0], 0],
    #                                        [0, 0, 0]]))
    F = F_cyl
    # Right Cauchy-Green tensor
    C = ufl.variable(F.T * F)

    # Finger strain tensor
    f = ufl.variable(ufl.inv(C))

    # Green strain tensor
    B = ufl.variable(F*F.T)

    # piola strain tensor
    c = ufl.variable(ufl.inv(B))

    # Lagrangian finite strain tensor
    E_GL = ufl.variable((C-I)/2)

    # Lagrangian finite strain tensor
    e_EA = ufl.variable((I - ufl.inv(B)) / 2)

    # Invariants of deformation tensors
    Ic = ufl.variable(ufl.tr(C))
    J = ufl.variable(ufl.det(F))

    def rule_of_mix(par):
        return

    def eps(v):
        e = ufl.sym(ufl.grad(v))
        return e

    def eps_th(vTemp):
        e = alpha*(vTemp-T)
        return e

    def as_3D_tensor(X):
        return ufl.as_tensor([[X[0], X[3], 0], [X[3], X[1], 0], [0, 0, X[2]]])

    # Define the elasticity model via a stored strain energy density function $\psi$, and create the expression for the first Piola-Kirchhoff stress:

    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2

    # Hyper-elasticity
    P = ufl.diff(psi, F)
    Pi = psi * ufl.dx
    # Push to Cauchy stress
    sigma = ufl.variable(ufl.inv(J) * P * F.T)

    gdim = domain.geometry.dim

    u_sol = fem.Expression(u, V.element.interpolation_points())
    u_sol_val = fem.Function(V)
    u_sol_val.name = "Deformation"
    u_sol_val.interpolate(u_sol)


    eps_el_expr = fem.Expression(ufl.as_tensor(E_GL), VT.element.interpolation_points())
    eps_el_val = fem.Function(VT)
    eps_el_val.name = "Elastic strain"
    eps_el_val.interpolate(eps_el_expr)

    sigma_expr = fem.Expression(ufl.as_tensor(sigma.expression()), VT.element.interpolation_points())
    sigma_val = fem.Function(VT)
    sigma_val.name = "Stress"
    sigma_val.interpolate(sigma_expr)

    test_expr = fem.Expression(ufl.as_tensor(sigma.expression()), VT.element.interpolation_points())
    test_val = fem.Function(VT)
    test_val.name = "Stress"
    test_val.interpolate(sigma_expr)

    def get_values_scalar(var):
        tmpexpr = fem.Expression(var, Vscal.element.interpolation_points())
        var_val = fem.Function(Vscal)
        var_val.interpolate(tmpexpr)
        print(var_val.x.array)
        return var_val

    def get_values_vector(var):
        tmpexpr = fem.Expression(var, V.element.interpolation_points())
        var_val = fem.Function(V)
        var_val.interpolate(tmpexpr)
        print(var_val.x.array.reshape(-1, gdim,))
        return var_val
    def get_values(var):
        tmpexpr = fem.Expression(var, VT.element.interpolation_points())
        var_val = fem.Function(VT)
        var_val.interpolate(tmpexpr)
        tmp = var_val.x.array.reshape(-1, gdim, gdim)
        print(tmp)
        print(np.shape(tmp))
        return var_val

    # ```{admonition} Comparison to linear elasticity
    # To illustrate the difference between linear and hyperelasticity, the following lines can be uncommented to solve the linear elasticity problem.
    # ```




    # Define the variational form with traction integral over all facets with value 2. We set the quadrature degree for the integrals to 4.

    metadata = {"quadrature_degree": 2}
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)
    #dr = ufl.Measure("dr", domain=domain, metadata=metadata)
    #print(f"Cylinder volume: {2 * np.pi * fem.assemble(fem.Constant(domain, 1.0) * x[0] * dx(domain=domain))}")

    # Define form F (we want to find u such that F(u) = 0)
    Form = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, BodF) * dx
    # Cylindrical form?
    #Form = ufl.inner(sigma(u), epsilon(v)) * ufl.dx - ufl.dot(f, v) * ufl.dx
    F_form = ufl.derivative(Pi, u, v)
    J_form = ufl.derivative(F_form, u)
    # Form = 2 * np.pi * ufl.inner(ufl.grad(v), P) * x[0] * dx - 2 * np.pi * ufl.inner(v, BodF) * x[0] * dx
    # As the varitional form is non-linear and written on residual form, we use the non-linear problem class from DOLFINx to set up required structures to use a Newton solver.
    problem = NonlinearProblem(F_form, u, bcs, J=J_form)
    problem = NonlinearProblem(Form, u, bcs)

    # and then create and customize the Newton solver

    solver = NewtonSolver(domain.comm, problem)

    # Set Newton solver options
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"

    # We create a function to plot the solution at each time step.

    # Compute magnitude of displacement to visualize in GIF
    magnitude = fem.Function(Vscal)
    us = fem.Expression(ufl.sqrt(sum([u[i] ** 2 for i in range(len(u))])), Vscal.element.interpolation_points())
    magnitude.interpolate(us)
    # -

    #stress = fem.functionspace(V, "P", 2)
    #stress.interpolate(stress_exp)
    # Finally, we solve the problem over several time steps, updating the z-component of the traction

    log.set_log_level(log.LogLevel.INFO)
    tval0 = 1e9
    for n in range(1, 10):
        BodF.value[0] = n * tval0
        num_its, converged = solver.solve(u)
        assert (converged)
        u.x.scatter_forward()
        print(f"Time step {n}, Number of iterations {num_its}, Load {BodF.value}")
        sigma_val.interpolate(sigma_expr)
        eps_el_val.interpolate(eps_el_expr)
        u_sol_val.interpolate(u_sol)
        magnitude.interpolate(us)
        test_val(test_expr)
        #print(stress.x.array)
        #print(fem.Expression(sig(eps(u)), u.x).eval(domain))
        #print(magnitude.x.array)
    #get_values(E_GL)
    #get_values(sigma)
    #get_values_vector(u)

    from pathlib import Path
    results_folder = Path("FeniCSx")
    results_folder.mkdir(exist_ok=True, parents=True)
    filename = results_folder / "beam_stress"
    results = [sigma_val, eps_el_val]
    with io.VTXWriter(domain.comm, filename.with_suffix(".bp"), results) as vtx:
        vtx.write(0.0)

    #with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "beam_stress.xdmf", "w") as xdmf:
    #    xdmf.write_mesh(domain)
    #    #xdmf.write_function(u_sol_val, "Displacement")
    #    xdmf.write_function(sigma_val, "Stress")
    #with io.XDMFFile(domain.comm, filename.with_suffix(".xdmf"), "w") as xdmf:
    #    xdmf.write_mesh(domain)
    #    xdmf.write_function(test_val)
    print("Data")
    #print(eps_el_val.x.array.reshape(-1, gdim, gdim)[2::2, :, :])
    #print(sigma_val.x.array.reshape(-1, gdim, gdim)[2::2, :, :])
    sigma_T = sigma_val.x.array.reshape(-1, gdim, gdim)
    elementnr = 4
    indx = np.array([0, 2, 1], dtype='int')
    for i in range(elementnr-1):
        indx = np.concatenate([indx, np.add(int(2*i), [4, 3])])
        print(indx)
    print(eps_el_val.x.array[0::9][indx])
    print(sigma_val.x.array[0::9][indx])
    print(sigma_val.x.array[4::9][indx])
    print(sigma_val.x.array[8::9][indx])
    print(sigma_T[indx])
    #print(fem.assemble(E_GL))
    #print(np.shape(sigma_val.x.array.reshape(-1, gdim, gdim)))
    #print(V.element.interpolation_points())
    geom = fem.Function(V)
    geom_exp = fem.Expression(x, V.element.interpolation_points())
    geom.interpolate(geom_exp)
    #print(geom.x.array.reshape(-1, gdim)[2::2, :])
    #sigma_val.x.array.reshape(-1, gdim, gdim)[2::2, :, :]
    #adjustxdmf()
    import meshio

    with meshio.xdmf.TimeSeriesReader("FCSX.xdmf") as reader:
        points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            datanames = point_data.keys()
            print(datanames)
def FNXTest_1D():
    # +
    from dolfinx import log, default_scalar_type
    from dolfinx.fem.petsc import NonlinearProblem
    from dolfinx.nls.petsc import NewtonSolver
    import numpy as np
    import ufl

    from mpi4py import MPI
    from dolfinx import fem, mesh, plot
    import basix

    # Creating calculation grid
    gmshmodel = gmsh1D_Really1D()
    createxdmf()

    gdim = 1
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)
    V = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    VT = fem.functionspace(domain, ("Lagrange", 2, (3, 3)))
    Vscal = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    x = SpatialCoordinate(domain)

    geom = fem.Function(V)
    geom_exp = fem.Expression(x, V.element.interpolation_points())
    geom.interpolate(geom_exp)
    elementnr = 4
    indx = np.array([0, 2, 1], dtype='int')
    for i in range(elementnr - 1):
        indx = np.concatenate([indx, np.add(int(2 * i), [4, 3])])
    print(geom.x.array[indx])


    # Elasticity parameters
    E = default_scalar_type(200.0e9)
    nu = default_scalar_type(0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
    alpha = default_scalar_type(2.2e-05)
    T0 = default_scalar_type(1073.15)
    TRef = default_scalar_type(1073.15)


    # -

    # We then create a function for supplying the boundary condition on the left side, which is fixed.

    u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)

    # To apply the boundary condition, we identity the dofs located on the facets marked by the `MeshTag`.
    def center(x):
        return np.isclose(x[0], 0)

    def circ(x):
        return np.isclose(x, 0.008)

    center_dofs = fem.locate_dofs_geometrical(V, center)
    #circ_dofs = fem.locate_dofs_geometrical(V, circ)
    #center_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
    # V.sub(0)
    #bcs = [fem.dirichletbc(u_bc, center_dofs, V), fem.dirichletbc(u_bc, circ_dofs, V)]
    bcs = [fem.dirichletbc(u_bc, center_dofs, V)]

    # Next, we define the body force on the reference configuration (`B`), and nominal (first Piola-Kirchhoff) traction (`T`).
    BodF = fem.Constant(domain, default_scalar_type((0,)))
    Trac = fem.Constant(domain, default_scalar_type(0,))

    # Define the test and solution functions on the space $V$

    u = fem.Function(V)
    u.name = "Displacement"
    v = ufl.TestFunction(V)

    Temp = fem.Function(V)
    Temp.name = "Temperature"
    dTemp = ufl.TestFunction(V)

    fM = fem.Function(Vscal)
    fM.name = "Martensite"
    dfM = ufl.TestFunction(Vscal)

    fB = fem.Function(Vscal)
    fB.name = "Bainate"
    dfB = ufl.TestFunction(Vscal)

    fP = fem.Function(Vscal)
    fP.name = "Pearlite"
    dfP = ufl.TestFunction(Vscal)

    fA = fem.Function(Vscal)
    fA.name = "Austenite"
    dfA = ufl.TestFunction(Vscal)




    # Define kinematic quantities used in the problem

    # Spatial dimension
    d = len(u)

    # Identity tensor
    I = ufl.variable(ufl.Identity(3))

    # Deformation gradient
    #F = ufl.variable(I + ufl.grad(u))
    gradxyz = ufl.grad(u)
    #F_cyl = ufl.variable(I + ufl.as_tensor([[gradxyz[0, 0], (gradxyz[0, 1] - u[1])/x[0], gradxyz[0, 2]],
    #                                        [gradxyz[1, 0], (gradxyz[1, 1] + u[0])/x[0], gradxyz[1, 2]/x[0]],
    #                                        [gradxyz[2, 0], gradxyz[2, 1]/x[0], gradxyz[2, 2]]]))
    F_cyl = ufl.variable(I + ufl.as_tensor([[gradxyz[0, 0], 0,          0],
                                            [0,             u[0]/x[0],  0],
                                            [0,             0,          -nu/(1-nu)*(gradxyz[0, 0] + u[0]/x[0])]]))
    #F_cyl = ufl.variable(I + ufl.as_tensor([[gradxyz[0, 0], (gradxyz[0, 1] - u[1]) / x[0], 0],
    #                                        [gradxyz[1, 0], (gradxyz[1, 1] + u[0]) / x[0], 0],
    #                                        [0, 0, 0]]))
    F = F_cyl
    # Right Cauchy-Green tensor
    C = ufl.variable(F.T * F)

    # Finger strain tensor
    f = ufl.variable(ufl.inv(C))

    # Green strain tensor
    B = ufl.variable(F*F.T)

    # piola strain tensor
    c = ufl.variable(ufl.inv(B))

    # Lagrangian finite strain tensor
    E_GL = ufl.variable((C-I)/2)

    # Lagrangian finite strain tensor
    e_EA = ufl.variable((I - ufl.inv(B)) / 2)

    # Invariants of deformation tensors
    Ic = ufl.variable(ufl.tr(C))
    J = ufl.variable(ufl.det(F))

    def rule_of_mix(par):
        return

    def eps(v):
        e = ufl.sym(ufl.grad(v))
        return e

    def eps_th(vTemp):
        e = alpha*(vTemp-T)
        return e

    def as_3D_tensor(X):
        return ufl.as_tensor([[X[0], X[3], 0], [X[3], X[1], 0], [0, 0, X[2]]])

    # Define the elasticity model via a stored strain energy density function $\psi$, and create the expression for the first Piola-Kirchhoff stress:

    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2

    # Hyper-elasticity
    P = ufl.diff(psi, F)

    # Push to Cauchy stress
    sigma = ufl.variable(ufl.inv(J) * P * F.T)

    gdim = domain.geometry.dim

    u_sol = fem.Expression(u, V.element.interpolation_points())
    u_sol_val = fem.Function(V)
    u_sol_val.name = "Deformation"
    u_sol_val.interpolate(u_sol)


    eps_el_expr = fem.Expression(ufl.as_tensor(E_GL), VT.element.interpolation_points())
    eps_el_val = fem.Function(VT)
    eps_el_val.name = "Elastic strain"
    eps_el_val.interpolate(eps_el_expr)

    sigma_expr = fem.Expression(ufl.as_tensor(sigma.expression()), VT.element.interpolation_points())
    sigma_val = fem.Function(VT)
    sigma_val.name = "Stress"
    sigma_val.interpolate(sigma_expr)

    test_expr = fem.Expression(ufl.as_tensor(sigma.expression()), VT.element.interpolation_points())
    test_val = fem.Function(VT)
    test_val.name = "Stress"
    test_val.interpolate(sigma_expr)

    def get_values_scalar(var):
        tmpexpr = fem.Expression(var, Vscal.element.interpolation_points())
        var_val = fem.Function(Vscal)
        var_val.interpolate(tmpexpr)
        print(var_val.x.array)
        return var_val

    def get_values_vector(var):
        tmpexpr = fem.Expression(var, V.element.interpolation_points())
        var_val = fem.Function(V)
        var_val.interpolate(tmpexpr)
        print(var_val.x.array.reshape(-1, gdim,))
        return var_val
    def get_values(var):
        tmpexpr = fem.Expression(var, VT.element.interpolation_points())
        var_val = fem.Function(VT)
        var_val.interpolate(tmpexpr)
        tmp = var_val.x.array.reshape(-1, gdim, gdim)
        print(tmp)
        print(np.shape(tmp))
        return var_val

    # ```{admonition} Comparison to linear elasticity
    # To illustrate the difference between linear and hyperelasticity, the following lines can be uncommented to solve the linear elasticity problem.
    # ```




    # Define the variational form with traction integral over all facets with value 2. We set the quadrature degree for the integrals to 4.

    metadata = {"quadrature_degree": 2}
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)
    #dr = ufl.Measure("dr", domain=domain, metadata=metadata)
    #print(f"Cylinder volume: {2 * np.pi * fem.assemble(fem.Constant(domain, 1.0) * x[0] * dx(domain=domain))}")

    # Define form F (we want to find u such that F(u) = 0)
    tmp = ufl.grad(v)
    print("Shape" + str(np.shape(tmp)))
    print(tmp)
    gradv = ufl.as_tensor([[tmp[0, 0], 0, 0], [0, 0, 0], [0, 0, 0]])
    gradv = ufl.as_tensor([[tmp[0, 0], 0,          0],
                           [0,             v[0]/x[0],  0],
                           [0,             0,          -nu/(1-nu)*(tmp[0, 0] + v[0]/x[0])]])
    testing = ufl.inner(gradv, P)
    testing2 = ufl.inner(v, BodF)
    print(testing)
    print(testing2)

    Form = ufl.inner(gradv, P) * x[0] * dx - ufl.inner(v, BodF) * x[0] * dx
    print(Form)

    # Cylindrical form?
    #Form = ufl.inner(sigma(u), epsilon(v)) * ufl.dx - ufl.dot(f, v) * ufl.dx
    #F_form = ufl.derivative(Pi, u, v)
    #J_form = ufl.derivative(F_form, u)
    # Form = 2 * np.pi * ufl.inner(ufl.grad(v), P) * x[0] * dx - 2 * np.pi * ufl.inner(v, BodF) * x[0] * dx
    # As the varitional form is non-linear and written on residual form, we use the non-linear problem class from DOLFINx to set up required structures to use a Newton solver.
    #problem = NonlinearProblem(F_form, u, bcs, J=J_form)
    problem = NonlinearProblem(Form, u, bcs)

    # and then create and customize the Newton solver

    solver = NewtonSolver(domain.comm, problem)

    # Set Newton solver options
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"

    # We create a function to plot the solution at each time step.

    # Compute magnitude of displacement to visualize in GIF
    #magnitude = fem.Function(Vscal)
    #us = fem.Expression(ufl.sqrt(sum([u[i] ** 2 for i in range(len(u))])), Vscal.element.interpolation_points())
    #magnitude.interpolate(us)
    # -

    #stress = fem.functionspace(V, "P", 2)
    #stress.interpolate(stress_exp)
    # Finally, we solve the problem over several time steps, updating the z-component of the traction

    log.set_log_level(log.LogLevel.INFO)
    tval0 = 1e9
    for n in range(1, 10):
        BodF.value[0] = n * tval0
        num_its, converged = solver.solve(u)
        assert (converged)
        u.x.scatter_forward()
        print(f"Time step {n}, Number of iterations {num_its}, Load {BodF.value}")
        sigma_val.interpolate(sigma_expr)
        eps_el_val.interpolate(eps_el_expr)
        u_sol_val.interpolate(u_sol)
        test_val(test_expr)
        #print(stress.x.array)
        #print(fem.Expression(sig(eps(u)), u.x).eval(domain))
        #print(magnitude.x.array)
    #get_values(E_GL)
    #get_values(sigma)
    #get_values_vector(u)
    from pathlib import Path
    results_folder = Path("FeniCSx")
    results_folder.mkdir(exist_ok=True, parents=True)
    filename = results_folder / "beam_stress"
    results = [sigma_val, eps_el_val]
    with io.VTXWriter(domain.comm, filename.with_suffix(".bp"), results) as vtx:
        vtx.write(0.0)

    #with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "beam_stress.xdmf", "w") as xdmf:
    #    xdmf.write_mesh(domain)
    #    #xdmf.write_function(u_sol_val, "Displacement")
    #    xdmf.write_function(sigma_val, "Stress")
    #with io.XDMFFile(domain.comm, filename.with_suffix(".xdmf"), "w") as xdmf:
    #    xdmf.write_mesh(domain)
    #    xdmf.write_function(test_val)
    print("Data")
    #print(eps_el_val.x.array.reshape(-1, gdim, gdim)[2::2, :, :])
    #print(sigma_val.x.array.reshape(-1, gdim, gdim)[2::2, :, :])
    sigma_T = sigma_val.x.array.reshape(-1, gdim, gdim)
    print(eps_el_val.x.array[0::9][indx])
    print(sigma_val.x.array[0::9][indx])
    print(sigma_val.x.array[4::9][indx])
    print(sigma_val.x.array[8::9][indx])
    print("Displacements")
    print(u_sol_val.x.array[indx])
    print(sigma_T[indx])

    geom = fem.Function(V)
    geom_exp = fem.Expression(x, V.element.interpolation_points())
    geom.interpolate(geom_exp)
    #print(geom.x.array.reshape(-1, gdim)[2::2, :])
    #sigma_val.x.array.reshape(-1, gdim, gdim)[2::2, :, :]
    #adjustxdmf()
    sigma_T = sigma_val.x.array.reshape(-1, gdim, gdim)


    test_V = fem.functionspace(domain, ("Lagrange", 1, (3, 3)))

    sigma_expr = fem.Expression(ufl.as_tensor(sigma.expression()), test_V.element.interpolation_points())
    testinterpolation = fem.Function(test_V)
    testinterpolation.name = "Stress"
    testinterpolation.interpolate(sigma_expr)

    with io.VTXWriter(domain.comm, filename.with_suffix(".bp"), [testinterpolation]) as vtx:
        vtx.write(0.0)

    with io.XDMFFile(domain.comm, "FCSX.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        #sigma_T.name = "Stress"
        #xdmf.write_function(u_sol_val)
        xdmf.write_function(testinterpolation, 0.0)

    import meshio
    try:
        with meshio.xdmf.TimeSeriesReader("FCSX.xdmf") as reader:
            reader.read_data(0)
            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
        print("RESULTS")
        print(t_list)
        print(cd_list)
        print(pd_list)
    except meshio._exceptions.ReadError:
        raise KeyError("No meshgrid for timeseries")

    return
    # Example


    with io.XDMFFile(domain.comm, filename.with_suffix(".xdmf"), "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(uh)


    import meshio


    try:
        with meshio.xdmf.TimeSeriesReader("FCSX.xdmf") as reader:
            points, cells = reader.read_points_cells()
            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
    except meshio._exceptions.ReadError:
        raise KeyError("No meshgrid for timeseries")
    t_list.append(1.0)
    indx_nodes = [0,2,4,6,8]
    with meshio.xdmf.TimeSeriesWriter("FCSX.xdmf") as writer:
        writer.write_points_cells(points, cells)
        writer.write_data(t=1.0, point_data={"sigmaXX": sigma_val.x.array[0::9][indx][indx_nodes]})

if __name__ == "__main__":
    FNXTest_1D()