import numpy as np
import ufl
import gmsh
from mpi4py import MPI
from dolfinx import fem, io, nls, default_scalar_type, log
import dolfinx.fem.petsc
from dolfinx.fem.petsc import LinearProblem
import dolfinx.nls.petsc
from dolfinx.mesh import create_box, CellType
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
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

def create_mesh(data):

    print('Meshing for FeniCSx specifics')

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
    gmsh.model.mesh.setOrder(1)
    # ----------------------
    print(*gmsh.logger.get(), sep="\n")

    gmsh.write("FeniCSx/FCSX_Mesh.msh")
    print("Created 1D mesh for FeniCSx")

    return gmsh.model

def create_xdmf():
    print('Creating temporary XDMF for FeniCSx specifics')
    filename = "FCSX.xdmf"
    mesh = meshio.read("FeniCSx/FCSX_Mesh.msh")
    with meshio.xdmf.TimeSeriesWriter(filename) as writer:
        writer.write_points_cells(points=mesh.points, cells={"line": mesh.get_cells_type("line")})
    print("Created temporary XDMF for FeniCSx")

def cyl_grad(v, x, nu):
    return ufl.as_tensor([[ufl.grad(v)[0, 0], 0, 0],
                           [0, v[0] / x[0], 0],
                           [0, 0, -nu / (1 - nu) * (ufl.grad(v)[0, 0] + v[0] / x[0])]])

def Run_1D():
    data = dict()
    data["Geometry"] = {"nodes": 5, "meshscaling": 1, "radius": 0.008}

    gmshmodel = create_mesh(data)
    create_xdmf()

    gdim = 1
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)

    # Creating functionspaces for calculation points
    # DG -> interpolation points
    # Lagrange -> Node points
    V = fem.functionspace(domain, ("Lagrange", 2, (1,)))
    VT = fem.functionspace(domain, ("Lagrange", 2, (3, 3)))

    x = SpatialCoordinate(domain)

    # Elasticity parameters
    E = default_scalar_type(200.0e9)
    nu = default_scalar_type(0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))

    # Temperature dependance parameters
    alpha = default_scalar_type(2.2e-05)
    T0 = default_scalar_type(1073.15)
    TRef = default_scalar_type(1073.15)

    # Fixed boundary condition
    u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)

    def center(dofs):
        return np.isclose(dofs[0], 0)

    def circ(dofs):
        return np.isclose(dofs, data['Geometry']['radius'])

    # Double check if V not having dofs at center creates error
    center_dofs = fem.locate_dofs_geometrical(V, center)
    bcs = [fem.dirichletbc(u_bc, center_dofs, V)]

    # Body forces acting
    BodF = fem.Constant(domain, default_scalar_type((0,)))

    # Define the test and solution functions on the space $V$
    u = fem.Function(V)
    u.name = "Displacement"
    v = ufl.TestFunction(V)

    Temp = fem.Function(V)
    Temp.name = "Temperature"
    dTemp = ufl.TestFunction(V)

    fM = fem.Function(V)
    fM.name = "Martensite"
    dfM = ufl.TestFunction(V)

    fB = fem.Function(V)
    fB.name = "Bainate"
    dfB = ufl.TestFunction(V)

    fP = fem.Function(V)
    fP.name = "Pearlite"
    dfP = ufl.TestFunction(V)

    fA = fem.Function(V)
    fA.name = "Austenite"
    dfA = ufl.TestFunction(V)

    # ----------------------------------------------------------------------------------------------------------------#
    # Defining kinematic quantities used in the problem

    I = ufl.variable(ufl.Identity(3))

    # Deformation gradient
    gradxyz = ufl.grad(u)
    F_cyl = ufl.variable(I + ufl.as_tensor([[gradxyz[0, 0], 0, 0],
                                            [0, u[0] / x[0], 0],
                                            [0, 0, -nu / (1 - nu) * (gradxyz[0, 0] + u[0] / x[0])]]))
    F = F_cyl

    # Right Cauchy-Green tensor
    C = ufl.variable(F.T * F)

    # Finger strain tensor
    f = ufl.variable(ufl.inv(C))

    # Green strain tensor
    B = ufl.variable(F * F.T)

    # piola strain tensor
    c = ufl.variable(ufl.inv(B))

    # Lagrangian finite strain tensor
    E_GL = ufl.variable((C - I) / 2)

    # Lagrangian finite strain tensor
    e_EA = ufl.variable((I - ufl.inv(B)) / 2)

    # Invariants of deformation tensors
    Ic = ufl.variable(ufl.tr(C))
    J = ufl.variable(ufl.det(F))

    # Define the elasticity model via a stored strain energy density function $\psi$, and create the expression for the first Piola-Kirchhoff stress:

    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2

    # Hyper-elasticity
    P = ufl.diff(psi, F)

    # Push to Cauchy stress
    sigma = ufl.inv(J) * P * F.T

    Vnodes = fem.functionspace(domain, ("CG", 1, (3, 3)))
    Vintpoints = fem.functionspace(domain, ("DG", 0, (3, 3)))

    # Elastic strain tensor
    eps_el_exp = fem.Expression(ufl.as_tensor(E_GL), Vnodes.element.interpolation_points())
    eps_el_int = fem.Function(Vnodes)
    eps_el_int.name = "Elastic strain"
    eps_el_int.interpolate(eps_el_exp)

    # Stress tensor
    sigma_exp = fem.Expression(ufl.as_tensor(sigma), Vnodes.element.interpolation_points())
    sigma_int = fem.Function(Vnodes)
    sigma_int.name = "Stress"
    sigma_int.interpolate(sigma_exp)

    # print(np.shape(sigma_int.x.array))
    # print([tuple(j) for j in sigma_int.x.array.reshape(-1, 9)[:, [0, 4, 8]]])
    # Double check
    metadata = {"quadrature_degree": 2}
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)

    # Define form F (we want to find u such that F(u) = 0)
    Form = ufl.inner(cyl_grad(v, x, nu), P) * x[0] * dx - ufl.inner(v, BodF) * x[0] * dx

    # As the varitional form is non-linear and written on residual form, we use the non-linear problem class from DOLFINx to set up required structures to use a Newton solver.
    problem = NonlinearProblem(Form, u, bcs)

    solver = NewtonSolver(domain.comm, problem)

    # Set Newton solver options
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"

    # Opening XDMF file for printing results
    with meshio.xdmf.TimeSeriesReader("FCSX.xdmf") as reader:
        points, cells = reader.read_points_cells()

    with meshio.xdmf.TimeSeriesWriter("FCSX.xdmf", data_format='HDF') as writer:
        writer.write_points_cells(points, cells)

        log.set_log_level(log.LogLevel.INFO)
        tval0 = 1e9
        for n in range(1, 10):
            BodF.value[0] = n * tval0
            num_its, converged = solver.solve(u)
            assert (converged)
            u.x.scatter_forward()
            print(f"Time step {n}, Number of iterations {num_its}, Load {BodF.value}")
            sigma_int.interpolate(sigma_exp)
            writer.write_data(float(n), point_data={"stress": sigma_int.x.array.reshape(-1, 9)})


if __name__ == "__main__":
    Run_1D()