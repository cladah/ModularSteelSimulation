
import numpy as np
import ufl
import gmsh
from mpi4py import MPI
from dolfinx import fem, io, nls
import dolfinx.fem.petsc
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
from dolfinx.io import XDMFFile, gmshio
import gmsh  # type: ignore
def gmsh1D():

    print('Remeshing 1D')

    data = dict()
    data["Geometry"] = {"nodes": 50, "meshscaling": 0.9, "radius": 0.008}

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
    print("Created 1D mesh for FenicsX")

    return gmsh.model

def create_mesh(comm: MPI.Comm, model: gmsh.model, name: str, filename: str, mode: str):
    """Create a DOLFINx from a Gmsh model and output to file.

    Args:
        comm: MPI communicator top create the mesh on.
        model: Gmsh model.
        name: Name (identifier) of the mesh to add.
        filename: XDMF filename.
        mode: XDMF file mode. "w" (write) or "a" (append).
    """
    msh, ct, ft = gmshio.model_to_mesh(model, comm, rank=0)
    msh.name = name
    ct.name = f"{msh.name}_cells"
    ft.name = f"{msh.name}_facets"
    with XDMFFile(msh.comm, filename, mode) as file:
        msh.topology.create_connectivity(2, 3)
        file.write_mesh(msh)
        file.write_meshtags(
            ct, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )
        file.write_meshtags(
            ft, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )


def FNX():
    gmshmodel = gmsh1D()
    comm = MPI.COMM_WORLD
    gdim = 1
    mode = "w"
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    msh, ct, ft = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)
    filename = "Testing.xdmf"
    msh.name = "1DmeshCyl"
    ct.name = f"{msh.name}_cells"
    ft.name = f"{msh.name}_facets"
    with XDMFFile(msh.comm, filename, mode) as file:
        msh.topology.create_connectivity(2, 3)
        file.write_mesh(msh)
        file.write_meshtags(
            ct, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )
        file.write_meshtags(
            ft, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )
    return

    tmpfile = io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf", "r")
    mesh = tmpfile.read_mesh()
    #print(mesh)
    #print(test)

    dim = mesh.topology.dim
    print(f"Mesh topology dimension d={dim}.")
    print(f"Mesh geometry dimension d={mesh.geometry.dim}.")
    degree = 2
    shape = (dim,)
    gdim = dim
    v_cg2 = element("Lagrange", mesh.topology.cell_name(), degree, shape=(mesh.topology.dim, ))
    V = fem.functionspace(mesh, v_cg2)

    # Adjusting boundary conditions
    Vx, _ = V.sub(0).collapse()
    Vy, _ = V.sub(1).collapse()
    #bottom_dofsy = fem.locate_dofs_topological((V.sub(1), Vy), gdim - 1, facets.find(1))
    #top_dofsx = fem.locate_dofs_topological((V.sub(0), Vx), gdim - 1, facets.find(2))



    #V = fem.functionspace(mesh, ("P", degree, shape))

    u = fem.Function(V, name="Displacement")
    fT = fem.Function(V, name="Temperature")
    ffa = fem.Function(V, name="Austenite")
    ffp = fem.Function(V, name="Pearlite")
    ffb = fem.Function(V, name="Bainite")
    fff = fem.Function(V, name="Ferrite")
    ffM = fem.Function(V, name="Martensite")
    cC = fem.Function(V, name="C_Carbon")
    cN = fem.Function(V, name="C_N")
    # -

    W0e = element("Lagrange", mesh.topology.cell_name(), degree, shape=(mesh.topology.dim, ))
    We = element("Lagrange", mesh.topology.cell_name(), degree, shape=(mesh.topology.dim, ))
    W = fem.functionspace(mesh, We)
    W0 = fem.functionspace(mesh, W0e)

    sig = fem.Function(W)
    sig_old = fem.Function(W)
    n_elas = fem.Function(W)
    beta = fem.Function(W0)
    p = fem.Function(W0, name="Cumulative_plastic_strain")
    dp = fem.Function(W0)
    u = fem.Function(V, name="Total_displacement")
    du = fem.Function(V, name="Iteration_correction")
    Du = fem.Function(V, name="Current_increment")


    # Next, we define the corresponding hyperelastic potential using UFL operators. We can easily obtain the UFL expression for the PK1 stress by differentiating the potential $\psi$ {eq}`psi-neo-Hookean` with respect to the deformation gradient $\bF$. We therefore declare it as a variable using `ufl.variable` and then compute $\bP = \dfrac{\partial \psi}{\partial \bF}$ using `ufl.diff`.

    # +
    # Identity tensor
    Id = Identity(dim)
    # Deformation gradient
    F = variable(Id + grad(u))

    # Right Cauchy-Green tensor
    C = F.T * F

    # Invariants of deformation tensors
    I1 = tr(C)
    J = det(F)

    # Shear modulus
    E = 1e4
    nu = 0.3
    mu = fem.Constant(mesh, E / 2 / (1 + nu))
    lmbda = fem.Constant(mesh, E * nu / (1 - 2 * nu) / (1 + nu))
    sig0 = fem.Constant(mesh, 250.0)  # yield strength in MPa
    Et = E / 100.0  # tangent modulus
    H = E * Et / (E - Et)  # hardening modulus

    # Stored strain energy density (compressible neo-Hookean model)
    psi = mu / 2 * (I1 - 3 - 2 * ln(J)) + lmbda / 2 * (J - 1) ** 2

    # PK1 stress = d_psi/d_F
    P = ufl.diff(psi, F)
    print(P)


    # -

    # Now, we set up the boundary conditions by first identifying the top and bottom dofs. We use Functions to provide the imposed displacement on both faces. For now, such functions are zero.

    # +
    def bottom(x):
        return np.isclose(x[2], 0.0)


    def top(x):
        return np.isclose(np.arctan(x[2]/x[1]), np.pi/6)


    bottom_dofs = fem.locate_dofs_geometrical(V, bottom)
    top_dofs = fem.locate_dofs_geometrical(V, top)

    u_bot = fem.Function(V)
    u_top = fem.Function(V)

    bcs = [fem.dirichletbc(u_bot, bottom_dofs), fem.dirichletbc(u_top, top_dofs)]
    # -

    # We will later update the value of the `u_top` function based on a UFL expression corresponding to the imposed rigid body rotation. This expression depends on a scalar value $\theta$ represented as a `Constant` object. The use of a `fem.Expression` results in JIT compilation of the code corresponding to the evaluation of this expression at specific points in the reference elements (here the interpolation points of $V$ i.e. the hexahedron vertices).

    x = SpatialCoordinate(mesh)
    theta = fem.Constant(mesh, 0.0)
    # Rot = as_matrix([[cos(theta), sin(theta), 0], [-sin(theta), cos(theta), 0], [0, 0, 1]])
    # Changed the rotation matrix to 2D
    Rot = as_matrix([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
    rotation_displ = dot(Rot, x) - x
    rot_expr = fem.Expression(rotation_displ, V.element.interpolation_points())

    # Now, we define the global non-linear potential energy. Note that since we have non-linear expressions, we specify to the measure `dx` the desired level of accuracy of the quadrature method. Otherwise, FEniCS may use overly conservative estimates of the required number of quadrature points.
    #
    # Next, we compute the corresponding non-linear residual using the `ufl.derivative` function which computes the directional derivative in the direction of the TestFunction `v`.
    # We also apply it to the residual itself to compute the corresponding consistent tangent bilinear form, usually called the Jacobian in the context of a Newton method. The latter is computed in the direction of the TrialFunction `du`.

    # +
    dx = ufl.Measure("dx", domain=mesh, metadata={"quadrature_degree": 4})
    E_pot = psi * dx

    v = TestFunction(V)
    du = TrialFunction(V)
    Residual = derivative(
        E_pot, u, v
    )  # This is equivalent to Residual = inner(P, grad(v))*dx
    Jacobian = derivative(Residual, u, du)
    # -

    # Finally, we set up a `NonlinearProblem` instance based on the corresponding residual and jacobian, unknown function and boundary conditions. The latter will also be attached to a nonlinear solver implementing a Newton method.

    # +
    problem = fem.petsc.NonlinearProblem(Residual, u, bcs)

    solver = nls.petsc.NewtonSolver(mesh.comm, problem)
    # Set Newton solver options
    solver.atol = 1e-4
    solver.rtol = 1e-4
    solver.convergence_criterion = "incremental"
    # -

    # We are now in position to write the load-stepping loop which simply updates the value of $\theta$. Since, `rot_expr` is symbolically linked to `theta`, this new value is automatically accounted for when interpolating the imposed top displacement from `rot_expr`.

    # +
    angle_max = 2 * np.pi
    Nsteps = 15

    out_file = "hyperelasticity.xdmf"
    with io.XDMFFile(mesh.comm, out_file, "w") as xdmf:
        xdmf.write_mesh(mesh)

    u.vector.set(0.0)
    for n, angle in enumerate(np.linspace(0, angle_max, Nsteps + 1)[1:]):
        theta.value = angle
        u_top.interpolate(rot_expr)

        num_its, converged = solver.solve(u)
        assert converged

        u.x.scatter_forward()  # updates ghost values for parallel computations

        print(
            f"Time step {n}, Number of iterations {num_its}, Angle {angle*180/np.pi:.0f} deg."
        )
        with io.XDMFFile(mesh.comm, out_file, "a") as xdmf:
            xdmf.write_function(u, n + 1)
# -

# ## Exercise
#
# Change the hyperelasticity model to another version of a compressible neo-Hookean model:
# \begin{equation*}
# \psi(\bF) = \dfrac{\mu}{2}(\overline{I}_1-3) + \left(\dfrac{\mu}{12}+\dfrac{\lambda}{8}\right)\left(J^2 + J^{-2}-2 \right)
# \end{equation*}
# where $\overline{I}_1 = \tr(\overline{\bC})$ with $\overline{\bC} = J^{-2/3}\bC$.
if __name__ == "__main__":
    FNX()