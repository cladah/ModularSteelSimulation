import os

from .ModuleStructure_file import CalcModule
from .Solvers.MeshSolvers import gmshsolver, pygmshsolver
import meshio


class Meshingmodule(CalcModule):
    def __init__(self):
        super().__init__("Meshing")

    def run(self):
        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")

            with meshio.xdmf.TimeSeriesReader("Datastream_Cache.xdmf") as reader:
                points, cells = reader.read_points_cells()

            # meshio.write("Datastream.xdmf",
            #              meshio.Mesh(points=points,
            #                          cells={"triangle": cells.get_cells_type("triangle")}))
            os.remove("Datastream.xdmf")
            os.remove("Datastream.h5")
            with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
                writer.write_points_cells(points=points, cells=cells)
            print("Meshing module done\n")
            return

        if self.program == "Gmsh":
            gmshsolver(self)
        elif self.program == "PyGmsh":
            pygmshsolver(self)
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
