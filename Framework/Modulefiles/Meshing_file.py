import os

from .ModuleStructure_file_new import CalcModule
from .Solvers.MeshSolvers import gmshsolver, pygmshsolver
import meshio


class Meshingmodule(CalcModule):
    def __init__(self, infile):
        infile = "Cachefiles/" + infile + ".json"
        super().__init__("Meshing", infile)

    def run(self):
        outstr = ["Meshing module\n",
                  "Radius: " + str(self.ginput["Geometry"]["radius"]),
                  "Number of nodes: " + str(self.ginput["Geometry"]["nodes"]),
                  "Mesh scaling factor: " + str(self.ginput["Geometry"]["meshscaling"]),
                  "\n---------------------------------------------------------------------\n"]
        for line in outstr:
            self.writeoutput(line)
            print(line)


        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")

            with meshio.xdmf.TimeSeriesReader("Datastream_Cache.xdmf") as reader:
                points, cells = reader.read_points_cells()

            # meshio.write("Datastream.xdmf",
            #              meshio.Mesh(points=points,
            #                          cells={"triangle": cells.get_cells_type("triangle")}))
            try:
                os.remove("Datastream.xdmf")
                os.remove("Datastream.h5")
            except FileNotFoundError:
                pass
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
