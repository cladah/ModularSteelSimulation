
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
def gmsh1D():

    print('Remeshing 1D')

    data = dict()
    data["Geometry"] = {"nodes": 3, "meshscaling": 0.9, "radius": 0.008}

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
    print("Created 1D mesh for FenicsX")

    return gmsh.model

def FNX():
    gmshmodel = gmsh1D()
    gdim = 3
    mode = "w"
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)
    filename = "Resultfiles/FNXMesh.xdmf"
    domain.name = "1DmeshCyl"
    ct.name = f"{domain.name}_cells"
    facets.name = f"{domain.name}_facets"
    with XDMFFile(domain.comm, filename, mode) as file:
        domain.topology.create_connectivity(0, 1)
        file.write_mesh(domain)
        file.write_meshtags(
            ct, domain.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{domain.name}']/Geometry"
        )
        file.write_meshtags(
            facets, domain.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{domain.name}']/Geometry"
        )


    dim = domain.topology.dim
    print(f"Mesh topology dimension d={dim}.")
    print(f"Mesh geometry dimension d={domain.geometry.dim}.")
    degree = 1
    shape = (dim,)
    gdim = dim
    #v_cg2 = ufl.element("Lagrange", msh.topology.cell_name(), degree, shape=(msh.topology.dim, ))
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))
    #boundary_facets = fem.locate_entities_boundary(msh, 1, lambda x: np.full(x.shape[1], True, dtype=bool))
    #center = fem.locate_dofs_topological((V, 0.0), gdim - 1, facets.find(1))
    #outer = fem.locate_dofs_topological((V, 0.008), gdim - 1, facets.find(1))
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

    W0e = element("Lagrange", domain.topology.cell_name(), degree, shape=(domain.topology.dim, ))
    We = element("Lagrange", domain.topology.cell_name(), degree, shape=(domain.topology.dim, ))
    W = fem.functionspace(domain, We)
    W0 = fem.functionspace(domain, W0e)
    sig = fem.Function(W)
    sig_old = fem.Function(W)
    n_elas = fem.Function(W)
    beta = fem.Function(W0)
    eps_pl = fem.Function(W0, name="Cumulative_plastic_strain")
    deps_pl = fem.Function(W0)
    u = fem.Function(V, name="Total_displacement")
    du = fem.Function(V, name="Iteration_correction")
    Du = fem.Function(V, name="Current_increment")


    # Next, we define the corresponding hyperelastic potential using UFL operators. We can easily obtain the UFL expression for the PK1 stress by differentiating the potential $\psi$ {eq}`psi-neo-Hookean` with respect to the deformation gradient $\bF$. We therefore declare it as a variable using `ufl.variable` and then compute $\bP = \dfrac{\partial \psi}{\partial \bF}$ using `ufl.diff`.
    B = fem.Constant(domain, default_scalar_type((1, 0, 0)))
    T = fem.Constant(domain, default_scalar_type((0, 0, 0)))

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

    # Shear modulus



    # Elasticity parameters
    E = default_scalar_type(1.0e4)
    nu = default_scalar_type(0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
    sig0 = fem.Constant(domain, 250.0)  # yield strength in MPa
    Et = E / 100.0  # tangent modulus
    H = E * Et / (E - Et)  # hardening modulus
    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    # Stress
    # Hyper-elasticity
    P = ufl.diff(psi, F)
    # ------------------------------#



    # -

    # Now, we set up the boundary conditions by first identifying the top and bottom dofs. We use Functions to provide the imposed displacement on both faces. For now, such functions are zero.
    #left_dofs = fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))

    def center(x):
        return np.isclose(x[0], 0.0)

    def circ(x):
        return np.isclose(x[0], 0.008)

    center_dofs = fem.locate_dofs_geometrical(V, center)
    circ_dofs = fem.locate_dofs_geometrical(V, circ)
    center_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
    circ_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
    circ_bc = np.array([0.001, 0, 0])
    bcs = [fem.dirichletbc(center_bc, center_dofs, V)] # , fem.dirichletbc(circ_bc, circ_dofs, V)


    x = SpatialCoordinate(domain)
    # theta = fem.Constant(domain, 0.0)
    # Rot = as_matrix([[cos(theta), sin(theta), 0], [-sin(theta), cos(theta), 0], [0, 0, 1]])
    # Changed the rotation matrix to 2D
    #Rot = as_matrix([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
    #rotation_displ = dot(Rot, x) - x
    #rot_expr = fem.Expression(rotation_displ, V.element.interpolation_points())

    # Now, we define the global non-linear potential energy. Note that since we have non-linear expressions, we specify to the measure `dx` the desired level of accuracy of the quadrature method. Otherwise, FEniCS may use overly conservative estimates of the required number of quadrature points.
    #
    # Next, we compute the corresponding non-linear residual using the `ufl.derivative` function which computes the directional derivative in the direction of the TestFunction `v`.
    # We also apply it to the residual itself to compute the corresponding consistent tangent bilinear form, usually called the Jacobian in the context of a Newton method. The latter is computed in the direction of the TrialFunction `du`.

    # +
    v = TestFunction(V)
    du = fem.Function(V)

    metadata = {"quadrature_degree": 4}
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)
    #F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx
    E_pot = psi * dx
    Form = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx


    Residual = derivative(
        E_pot, u, v
    )  # This is equivalent to Residual = inner(P, grad(v))*dx
    #Jacobian = derivative(Residual, u, du)

    # Finally, we set up a `NonlinearProblem` instance based on the corresponding residual and jacobian, unknown function and boundary conditions. The latter will also be attached to a nonlinear solver implementing a Newton method.

    problem = fem.petsc.NonlinearProblem(Form, u, bcs)


    #a = ufl.dot(ufl.grad(du), ufl.grad(v)) * dx
    #L = p * v * dx
    #problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    #uh = problem.solve()

    solver = nls.petsc.NewtonSolver(domain.comm, problem)

    # Set Newton solver options
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"

    # -
    num_its, converged = solver.solve(u)
    print(u.eval(x, ct))
    assert converged
    u.x.scatter_forward()
    print(x)
    return

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

def FNXTest():
    # +
    from dolfinx import log, default_scalar_type
    from dolfinx.fem.petsc import NonlinearProblem
    from dolfinx.nls.petsc import NewtonSolver
    import numpy as np
    import ufl

    from mpi4py import MPI
    from dolfinx import fem, mesh, plot

    # Creating calculation grid
    gmshmodel = gmsh1D()
    gdim = 3
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)
    V = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))

    # We create two python functions for determining the facets to apply boundary conditions to

    fdim = domain.topology.dim - 1
    # -

    # We then create a function for supplying the boundary condition on the left side, which is fixed.

    u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)

    # To apply the boundary condition, we identity the dofs located on the facets marked by the `MeshTag`.
    def center(x):
        return np.isclose(x[0], 0)

    def circ(x):
        return np.isclose(x[0], 0.008)

    center_dofs = fem.locate_dofs_geometrical(V, center)
    center_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
    V.sub(0)
    bcs = [fem.dirichletbc(u_bc, center_dofs, V)]
    #bcs = [fem.dirichletbc(u_bc, center_dofs, V)]

    # Next, we define the body force on the reference configuration (`B`), and nominal (first Piola-Kirchhoff) traction (`T`).
    B = fem.Constant(domain, default_scalar_type((0, 0, 0)))
    T = fem.Constant(domain, default_scalar_type((0, 0, 0)))

    # Define the test and solution functions on the space $V$

    v = ufl.TestFunction(V)
    u = fem.Function(V)

    Temp = fem.Function(V)
    Temp.name = "Temperature"
    dTemp = ufl.TestFunction(V)

    sig = fem.Function(V)



    # Define kinematic quantities used in the problem

    # +
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
    # -

    # Define the elasticity model via a stored strain energy density function $\psi$, and create the expression for the first Piola-Kirchhoff stress:

    # Elasticity parameters
    E = default_scalar_type(1.0e4)
    nu = default_scalar_type(0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    # Stress
    # Hyper-elasticity
    P = ufl.diff(psi, F)

    # ```{admonition} Comparison to linear elasticity
    # To illustrate the difference between linear and hyperelasticity, the following lines can be uncommented to solve the linear elasticity problem.
    # ```

    def sig(eps_el):
        return lmbda * ufl.tr(eps_el) * ufl.Identity(3) + 2 * mu * eps_el

    def eps(v):
        e = ufl.sym(ufl.grad(v))
        return e


    # +
    # P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I
    # -

    # Define the variational form with traction integral over all facets with value 2. We set the quadrature degree for the integrals to 4.

    metadata = {"quadrature_degree": 4}
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)

    #f = fem.Constant(domain, default_scalar_type(0))
    #a = Temp * dTemp * dx + dt * ufl.dot(ufl.grad(Temp), ufl.grad(dTemp)) * dx
    #L = (u_n + dt * f) * v * dx
    #bilinear_form = fem.form(a)
    #linear_form = fem.form(L)

    # Define form F (we want to find u such that F(u) = 0)
    F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx

    # As the varitional form is non-linear and written on residual form, we use the non-linear problem class from DOLFINx to set up required structures to use a Newton solver.

    problem = NonlinearProblem(F, u, bcs)

    # and then create and customize the Newton solver

    # +
    solver = NewtonSolver(domain.comm, problem)

    # Set Newton solver options
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"

    # We create a function to plot the solution at each time step.

    # Compute magnitude of displacement to visualize in GIF
    Vs = fem.functionspace(domain, ("Lagrange", 2))
    magnitude = fem.Function(Vs)
    us = fem.Expression(ufl.sqrt(sum([u[i] ** 2 for i in range(len(u))])), Vs.element.interpolation_points())
    magnitude.interpolate(us)
    # -

    # Finally, we solve the problem over several time steps, updating the z-component of the traction

    log.set_log_level(log.LogLevel.INFO)
    tval0 = 1
    for n in range(1, 10):
        B.value[0] = n * tval0
        num_its, converged = solver.solve(u)
        assert (converged)
        u.x.scatter_forward()
        print(f"Time step {n}, Number of iterations {num_its}, Load {T.value}")
        print(u.x.array[0::6])

        print(fem.Expression(sig(eps(u)), u.x).eval(domain))
        #print(magnitude.x.array)

    # <img src="./deformation.gif" alt="gif" class="bg-primary mb-1" width="800px">
def Solverloop(Du, sig):
    while nRes / nRes0 > tol and niter < Nitermax:
        # solve for the displacement correction
        tangent_problem.assemble_lhs()
        tangent_problem.solve_system()

        # update the displacement increment with the current correction
        Du.vector.axpy(1, du.vector)  # Du = Du + 1*du
        Du.x.scatter_forward()

        # interpolate the new stresses and internal state variables
        interpolate_quadrature(sig_, sig)
        interpolate_quadrature(n_elas_, n_elas)
        interpolate_quadrature(beta_, beta)

        # compute the new residual
        tangent_problem.assemble_rhs()
        nRes = tangent_problem._b.norm()

        niter += 1
def FNXTest_2():
    # ---
    # jupyter:
    #   jupytext:
    #     formats: ipynb,py:light
    #     text_representation:
    #       extension: .py
    #       format_name: light
    #       format_version: '1.5'
    #       jupytext_version: 1.16.4
    #   kernelspec:
    #     display_name: Python 3 (ipykernel)
    #     language: python
    #     name: python3
    # ---

    # # Hyperelasticity
    # Author: Jørgen S. Dokken and Garth N. Wells
    #
    # This section shows how to solve the hyperelasticity problem for deformation of a beam.
    #
    # We will also show how to create a constant boundary condition for a vector function space.
    #
    # We start by importing DOLFINx and some additional dependencies.
    # Then, we create a slender cantilever consisting of hexahedral elements and create the function space `V` for our unknown.

    # +
    from dolfinx import log, default_scalar_type
    from dolfinx.fem.petsc import NonlinearProblem
    from dolfinx.nls.petsc import NewtonSolver
    import numpy as np
    import ufl

    from mpi4py import MPI
    from dolfinx import fem, mesh, plot
    L = 20.0
    domain = mesh.create_box(MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [L, 1, 1]], [2, 1, 1], mesh.CellType.hexahedron)
    gmshmodel = gmsh1D()
    gdim = 3
    mode = "w"
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)

    V = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))

    # -

    # We create two python functions for determining the facets to apply boundary conditions to

    # +
    def left(x):
        return np.isclose(x[0], 0)

    def right(x):
        return np.isclose(x[0], L)

    fdim = domain.topology.dim - 1
    left_facets = mesh.locate_entities_boundary(domain, fdim, left)
    right_facets = mesh.locate_entities_boundary(domain, fdim, right)
    # -

    # Next, we create a  marker based on these two functions

    # Concatenate and sort the arrays based on facet indices. Left facets marked with 1, right facets with two
    marked_facets = np.hstack([left_facets, right_facets])
    marked_values = np.hstack([np.full_like(left_facets, 1), np.full_like(right_facets, 2)])
    sorted_facets = np.argsort(marked_facets)
    facet_tag = mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

    # We then create a function for supplying the boundary condition on the left side, which is fixed.

    u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)

    # To apply the boundary condition, we identity the dofs located on the facets marked by the `MeshTag`.

    left_dofs = fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))
    bcs = [fem.dirichletbc(u_bc, left_dofs, V)]

    # Next, we define the body force on the reference configuration (`B`), and nominal (first Piola-Kirchhoff) traction (`T`).

    B = fem.Constant(domain, default_scalar_type((1, 0, 0)))
    T = fem.Constant(domain, default_scalar_type((0, 0, 0)))

    # Define the test and solution functions on the space $V$

    v = ufl.TestFunction(V)
    u = fem.Function(V)

    # Define kinematic quantities used in the problem

    # +
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
    # -

    # Define the elasticity model via a stored strain energy density function $\psi$, and create the expression for the first Piola-Kirchhoff stress:

    # Elasticity parameters
    E = default_scalar_type(1.0e4)
    nu = default_scalar_type(0.3)
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
    # Stress
    # Hyper-elasticity
    P = ufl.diff(psi, F)

    # ```{admonition} Comparison to linear elasticity
    # To illustrate the difference between linear and hyperelasticity, the following lines can be uncommented to solve the linear elasticity problem.
    # ```

    # +
    # P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I
    # -

    # Define the variational form with traction integral over all facets with value 2. We set the quadrature degree for the integrals to 4.

    metadata = {"quadrature_degree": 4}
    ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
    dx = ufl.Measure("dx", domain=domain, metadata=metadata)
    # Define form F (we want to find u such that F(u) = 0)
    F = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)

    # As the varitional form is non-linear and written on residual form, we use the non-linear problem class from DOLFINx to set up required structures to use a Newton solver.

    problem = NonlinearProblem(F, u, bcs)

    # and then create and customize the Newton solver

    # +
    solver = NewtonSolver(domain.comm, problem)

    # Set Newton solver options
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"

    # -

    # We create a function to plot the solution at each time step.

    # +
    # Compute magnitude of displacement to visualize in GIF
    Vs = fem.functionspace(domain, ("Lagrange", 2))
    magnitude = fem.Function(Vs)
    us = fem.Expression(ufl.sqrt(sum([u[i] ** 2 for i in range(len(u))])), Vs.element.interpolation_points())
    magnitude.interpolate(us)
    # -

    # Finally, we solve the problem over several time steps, updating the z-component of the traction

    log.set_log_level(log.LogLevel.INFO)
    tval0 = -1.5
    for n in range(1, 10):
        T.value[2] = n * tval0
        num_its, converged = solver.solve(u)
        assert (converged)
        u.x.scatter_forward()
        print(f"Time step {n}, Number of iterations {num_its}, Load {T.value}")
        magnitude.interpolate(us)
        # print(u.x.array)
        print(magnitude.x.array)

if __name__ == "__main__":
    FNXTest()