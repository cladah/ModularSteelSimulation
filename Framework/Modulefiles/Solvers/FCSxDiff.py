from mpi4py import MPI
import numpy as np
import sys
import os
from dolfinx import mesh as dmesh, fem, io, nls
import ufl
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
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from petsc4py import PETSc
from dolfinx.io import XDMFFile, gmshio
import adios4dolfinx
import dolfinx
from dolfinx.fem import Function, FunctionSpace, Constant, dirichletbc, form, locate_dofs_topological, assemble_scalar, assemble_vector
from Framework.HelpFile import readresultfile, read_modinput, read_geninput
import meshio

from Framework.Datastream_file import createdatastream, readdatastream

def FCSxDiffSolver(ginput, minput, mobility, compgrid):

    print("Running FeniCSx diffusion solver")
    width = ginput["Geometry"]["width"]
    height = ginput["Geometry"]["height"]


    boosts = minput["BoostNr"]
    boost_t = minput["BoostTime"]
    diff_t = minput["DiffTime"]
    kb = 1.380649E-23
    T = minput["Temp"][0]
    D = kb*T*mobility

    cwd = os.getcwd()
    ginput = read_geninput()
    mesh_filename = cwd + "/Resultfiles/" + ginput["Datastream"]["Savedirect"]
    comp = readdatastream("Composition/C")
    #domain = gmshio.read_from_msh("Resultfiles/Mesh.xdmf", MPI.COMM_WORLD, gdim=2)
    with io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf", "r") as infile:
        domain = infile.read_mesh(name="Grid")
    V = fem.functionspace(domain, ("Lagrange", 2))

    #x = SpatialCoordinate(domain)
    c_bound = fem.Constant(domain, 1.0)


    tdim = domain.topology.dim
    fdim = tdim - 1

    def bottom_bound(x):
        return np.isclose(x[1], 0.0)

    def top_bound(x):
        return np.isclose(x[1], height)

    top_facets = dolfinx.mesh.locate_entities_boundary(domain, fdim, top_bound)
    bc_top = fem.dirichletbc(c_bound, fem.locate_dofs_topological(V, fdim, top_facets), V)
    bottom_facets = dolfinx.mesh.locate_entities_boundary(domain, fdim, bottom_bound)
    bc_bottom = fem.dirichletbc(c_bound, fem.locate_dofs_topological(V, fdim, bottom_facets), V)

    bcs = [bc_top, bc_bottom]

    c_n1 = fem.Function(V)
    c_n = fem.Function(V)

    #print(np.shape(c_n.x.array[:]))
    #print(np.shape(compgrid["C"]))
    c_n.x.array[:] = compgrid["C"]
    print(c_n.x.array[:])
    #exit()
    v = ufl.TestFunction(V)
    dc = ufl.TrialFunction(V)
    dt = 0.01


    F = (c_n1 - c_n) / dt * v * ufl.dx + D * ufl.dot(ufl.grad(c_n1), ufl.grad(v)) * ufl.dx

    problem = fem.petsc.NonlinearProblem(F, c_n1, bcs)
    solver = nls.petsc.NewtonSolver(MPI.COMM_WORLD, problem)
    solver.rtol = 1e-8

    t = 0.0
    t_end = 1.0
    while t < t_end + 1e-8:
        t += dt

        n_iter, converged = solver.solve(c_n1)
        c_n1.x.scatter_forward()

        # Write solution to file

        # Update for next step
        c_n.x.array[:] = c_n1.x.array[:]

    compgrid["C"] = c_n.x.array[:]
    return compgrid