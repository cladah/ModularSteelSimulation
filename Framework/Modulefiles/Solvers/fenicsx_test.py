import numpy as np
import ufl

from mpi4py import MPI
from dolfinx import common, fem, io, nls, plot
import dolfinx.fem.petsc
import dolfinx.nls.petsc
from dolfinx.mesh import create_box, CellType
from petsc4py import PETSc
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

import dolfinx_mpc.utils
from dolfinx_mpc import LinearProblem, MultiPointConstraint
import json
# ---------------------------------------------------------#
# Reading datastream
# ---------------------------------------------------------#
print("Reading data from Datastream")
mesh = io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf", "r")
domain = mesh.read_mesh()

# ---------------------------------------------------------#
# Reading datastream
# ---------------------------------------------------------#
f = open('Cachefiles/input.json', 'r')
data = json.load(f)
f.close()

# ---------------------------------------------------------#
# Setting up function spaces
# ---------------------------------------------------------#

dim = domain.topology.dim
print(f"Mesh topology dimension d={dim}.")
print(f"Mesh geometry dimension d={domain.geometry.dim}.")
el_2 = element("Lagrange", domain.topology.cell_name(), 2, shape=(domain.topology.dim, ))
V = fem.functionspace(domain, el_2)

# ---------------------------------------------------------#
# Setting up functions
# ---------------------------------------------------------#

x = SpatialCoordinate(domain)
n = ufl.FacetNormal(domain)
u = fem.Function(V, name="displacement")
T = fem.Function(V, name="temperature")
ffa = fem.Function(V, name="Austenite")
ffp = fem.Function(V, name="Pearlite")
ffb = fem.Function(V, name="Bainite")
fff = fem.Function(V, name="Ferrite")
ffM = fem.Function(V, name="Martensite")
cC = fem.Function(V, name="C_C")
cN = fem.Function(V, name="C_N")
eps_pl = fem.Function(V, name="eps_pl")

# ---------------------------------------------------------#
# Setting up trail functions
# ---------------------------------------------------------#

du = TrialFunction(V)
dv = TestFunction(V)

# ---------------------------------------------------------#
# Setup continuum model
# ---------------------------------------------------------#

# Identity tensor
Id = Identity(dim)
# Deformation gradient
F = variable(Id + grad(u))
# Right Cauchy-Green tensor
C = F.T * F
# Left Cauchy-Green tensor
b = F * F.T

# Invariants of deformation tensors
I1 = tr(C)
I2 = (tr(C)*tr(C)-tr(C*C))/2
J = det(F)

# Material properties
E = data["Material"]["Austenite"]["E"]["E"][0]
h = data["Material"]["Austenite"]["h"]["h"][0]
""" Insert rule of mixture """ """ Some of the material properties have dependence on T """
nu = 0.3
mu = fem.Constant(domain, E / 2 / (1 + nu))
lmbda = fem.Constant(domain, E * nu / (1 - 2 * nu) / (1 + nu))

# ---------------------------------------------------------#
# Energy density
# ---------------------------------------------------------#

# Elastic
psi_e = mu / 2 * (I1 - 3 - 2 * ln(J)) + lmbda / 2 * (J - 1) ** 2
# Plastic

# Thermal

# Phase transformation

psi = psi_e

# ---------------------------------------------------------#
# Boundary conditions
# ---------------------------------------------------------#

def x_axis(x):
    return np.isclose(x[1], 0.0)


def wedge(x):
    return np.isclose(np.arctan(x[1]/x[0]), np.pi/6)

def circ(x):
    return np.isclose(np.sqrt(x[1]**2+x[0]**2), 0.01)

bottom_dofs = fem.locate_dofs_geometrical(V, x_axis)
wedge_dofs = fem.locate_dofs_geometrical(V, wedge)

n_diagonal = ufl.as_vector((np.cos(np.pi/6+np.pi/2), np.sin(np.pi/6+np.pi/2)))




#rotation_displ = x
#rol_expr = fem.Expression(rotation_displ, wedge_dofs)
#roller = fem.Expression(rol_expr, wedge_dofs)
#boundary_facets = fem.locate_entities_boundary(domain, domain.topology.dim - 1, x_axis)
#boundary_dofs_x = fem.locate_dofs_topological(V.sub(0), domain.topology.dim - 1, boundary_facets)

roll_bc = fem.Function(V)
with roll_bc.vector.localForm() as loc:
    loc.set(0.0)

def wedge_disp(x):
    x[1] = np.tan(np.pi/6)*x[0]
    return x
u_L = fem.Expression(x[1], V.element.interpolation_points())
#wedge_expr = fem.Expression(wedge_disp, V.element.interpolation_points())
#wedge_dofs.interpolate(wedge_expr)

bcs = [fem.dirichletbc(roll_bc, bottom_dofs), fem.dirichletbc(roll_bc, wedge_dofs)]

# Traction force
Trac = fem.Constant(domain, 500.)
# Body force
B = fem.Constant(domain, (0., 0.))


n = dolfinx_mpc.utils.create_normal_approximation(V, domain, 1)

with common.Timer("~Stokes: Create slip constraint"):
    mpc = MultiPointConstraint(V)
    mpc.create_slip_constraint(V, (domain, 1), n, bcs=bcs)
mpc.finalize()


# ---------------------------------------------------------#
# Variational problem
# ---------------------------------------------------------#

P = ufl.diff(psi, F)
dx = ufl.Measure("dx", domain=domain, metadata={"quadrature_degree": 4})
ds = ufl.Measure("ds", domain=domain)

# potential energy of the system
E_pot = psi * dx - dot(B, u) * dx - dot(Trac*n, u)*ds
# Minimizing of potential energy
Residual = derivative(E_pot, u, dv)  # This is equivalent to Residual = inner(P, grad(v))*dx
Jacobian = derivative(Residual, u, du)

# ---------------------------------------------------------#
# Setting up problem to solve
# ---------------------------------------------------------#

problem = fem.petsc.NonlinearProblem(Residual, u, bcs)

solver = nls.petsc.NewtonSolver(domain.comm, problem)
solver.atol = 1e-4
solver.rtol = 1e-4
solver.convergence_criterion = "incremental"

num_its, converged = solver.solve(u)
u.x.scatter_forward()

print(u.x.array)
print(np.min(u.x.array))
print(np.max(u.x.array))
