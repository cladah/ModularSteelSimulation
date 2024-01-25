from .ModuleStructure_file import CalcModule
from .Solvers.MeshSolvers import gmshsolver, pygmshsolver
import meshio


class Meshingmodule(CalcModule):
    def __init__(self):
        super().__init__("Meshing")

    def run(self):
        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")

            meshdata = meshio.read("Cachefiles/Datastream.xdmf")
            meshio.write("Resultfiles/Datastream.xdmf",
                         meshio.Mesh(points=meshdata.points,
                                     cells={"triangle": meshdata.get_cells_type("triangle")}))

            print("Meshing module done\n")
            return

        if self.program == "Gmsh":
            gmshsolver(self)
        elif self.program == "PyGmsh":
            pygmshsolver(self)
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
