import os

from .ModuleStructure_file_new import CalcModule
from .Solvers.MeshSolvers import gmshsolver, pygmshsolver
import meshio
from Framework.Datastream_file import adjustdatastream, readdatastream
import numpy as np

class Meshingmodule(CalcModule):
    def __init__(self, infile):
        infile = "Inputs/" + infile + ".json"
        super().__init__("Meshing", infile)

    def run(self):
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Meshing module: " + self.inputfile + "\n",
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
        print("Setting starting point for calculation")
        nodes = readdatastream("nodes")
        tmpelvalues = np.full(len(nodes), 0)
        for element in self.ginput["Material"]["Composition"].keys():
            tmpelvalues.fill(self.ginput["Material"]["Composition"][element])
            adjustdatastream({"Composition/" + element: tmpelvalues}, "nodes")


