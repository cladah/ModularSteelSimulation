import numpy as np
import ufl

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
from basix.ufl import element

# ---------------------------------------------------------#
# Reading datastream
# ---------------------------------------------------------#
print("Reading data from Datastream")
mesh = io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf", "r")
domain = mesh.read_mesh()

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
u = fem.Function(V, name="displacement")
T = fem.Function(V, name="temperature")

# ---------------------------------------------------------#
# Setting up trailfunctions
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
C = F * F.T

# Invariants of deformation tensors
I1 = tr(C)
I2 = (tr(C)*tr(C)-tr(C*C))/2
J = det(F)

# Material properties
E = 1e4
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
top_dofs = fem.locate_dofs_geometrical(V, wedge)

u_bot = fem.Function(V)
u_top = fem.Function(V)

bcs = [fem.dirichletbc(u_bot, bottom_dofs), fem.dirichletbc(u_top, top_dofs)]

press = 1.

# ---------------------------------------------------------#
# Variational problem
# ---------------------------------------------------------#

P = ufl.diff(psi, F)
dx = ufl.Measure("dx", domain=domain, metadata={"quadrature_degree": 4})
ds = ufl.Measure("ds", domain=domain, metadata={"quadrature_degree": 4})
# potential energy of the system
E_pot = psi * dx - press*ds
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