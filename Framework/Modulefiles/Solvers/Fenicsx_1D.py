
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

    x = SpatialCoordinate(domain)
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
    import basix

    # Creating calculation grid
    gmshmodel = gmsh1D()
    gdim = 3
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    domain, ct, facets = gmshio.model_to_mesh(gmshmodel, mesh_comm, gmsh_model_rank, gdim=gdim)
    V = fem.functionspace(domain, ("Lagrange", 2, (gdim,)))
    VT = fem.functionspace(domain, ("Lagrange", 1, (gdim, gdim)))
    Vscal = fem.functionspace(domain, ("Lagrange", 1, (1,)))
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
    bcs = [fem.dirichletbc(u_bc, center_dofs, V), fem.dirichletbc(u_bc, circ_dofs, V)]
    #bcs = [fem.dirichletbc(u_bc, center_dofs, V)]

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
                                            [gradxyz[1, 0], (gradxyz[1, 1] - u[0])/x[0], gradxyz[1, 2]/x[0]],
                                            [gradxyz[2, 0], gradxyz[2, 1]/x[0], gradxyz[2, 2]]]))
    F = F
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


    eps_el_expr = fem.Expression(eps(u), VT.element.interpolation_points())
    eps_el_val = fem.Function(VT)
    eps_el_val.name = "Elastic strain"
    eps_el_val.interpolate(eps_el_expr)

    sigma_expr = fem.Expression(sigma.expression(), VT.element.interpolation_points())
    sigma_val = fem.Function(VT)
    sigma_val.name = "Stress"
    sigma_val.interpolate(sigma_expr)

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

    # Define form F (we want to find u such that F(u) = 0)
    Form = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, BodF) * dx

    # As the varitional form is non-linear and written on residual form, we use the non-linear problem class from DOLFINx to set up required structures to use a Newton solver.

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
        #print(stress.x.array)
        #print(fem.Expression(sig(eps(u)), u.x).eval(domain))
        #print(magnitude.x.array)
    get_values(E_GL)
    #get_values(sigma)
    get_values_vector(u)

    from pathlib import Path
    results_folder = Path("FeniCSx")
    results_folder.mkdir(exist_ok=True, parents=True)
    filename = results_folder / "beam_stress"
    results = [sigma_val, eps_el_val]
    with io.VTXWriter(domain.comm, filename.with_suffix(".bp"), results) as vtx:
        vtx.write(0.0)

    # from adios2 import FileReader
    # from adios2 import Stream
    # import adios4dolfinx as adx
    return
    fullfile = filename.with_suffix(".bp")
    import adios2

    print(fullfile)
    with adios2.open(fullfile, "r") as fr:
        # Loop through each variable in the file
        for variable_name in fr.available_variables():
            # Inquire and read the variable
            data = fr.read(variable_name)
            print(f"Variable: {variable_name}")
            print(f"Data: {data}\n")
    return


    #with dolfinx.io.XDMFFile(MPI.COMM_WORLD, "beam_stress.xdmf", "w") as xdmf:
    #    xdmf.write_mesh(domain)
    #    #xdmf.write_function(u_sol_val, "Displacement")
    #    xdmf.write_function(sigma_val, "Stress")
    #with io.XDMFFile(domain.comm, filename.with_suffix(".xdmf"), "w") as xdmf:
    #    xdmf.write_mesh(domain)
    #    xdmf.write_function(u_sol_val)
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