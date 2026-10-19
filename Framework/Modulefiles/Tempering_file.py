from .ModuleStructure_file_new import CalcModule
import numpy as np
from numpy import interp
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream, getaxisvalues
from .Solvers.ThermocalcSolver import TCTemperingSolver

class Diffusionmodule(CalcModule):

    def __init__(self, infile, modulenr=0):
        infile = "Inputs/" + infile + ".json"
        super().__init__("Tempering", infile, modulenr)

    def run(self):
        Carbtime = sum(self.minput["BoostTime"] + self.minput["DiffTime"]) / 3600
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Tempering module: " + self.inputfile + "\n",
                  "\n---------------------------------------------------------------------\n"]

        for line in outstr:
            self.writeoutput(line)
            print(line)

        if not self.check_runcondition():
            print("NOT IMPLEMENTED")
            print("Using precalculated " + str(self.module) + " simulation")
            return

        if self.program == "TC":
            self.updateprogress(0.1)

            print('Running tempering module with ThermoCalc Prisma')

            self.updateprogress(0.2)

            composition = dict()
            for element in self.ginput["Material"]["Composition"].keys():
                composition[element] = getaxisvalues("Composition_" + element, -1)


            # Adding geometry to composition?
            composition = TCTemperingSolver(self.ginput, self.minput, composition)

            self.updateprogress(0.9)
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

        """
            Interpolating the compositional values to nodal points 
        """
        import pandas as pd
        df = pd.DataFrame()
        df["x"] = composition[0]
        df["C"] = composition[1]["C"]
        df.to_csv("CompositionData.csv")
        print("Interpolating result to nodal points")
        if self.ginput["Geometry"]["Type"] in ["4PointBend", "3PointBend"]:
            fourpointbend_interp(self, composition)
        elif self.ginput["Geometry"]["Type"] == "Cylinder":
            cylinder_interp(self, composition)
        else:
            raise KeyError("Type of geometry not supported for compositional "
                           "interpolation in diffusion modell")

        self.updateprogress(1.0)
        print("Carburization module done\n")

def fourpointbend_interp(parent, composition):
    print("interpolation for 4PB!")
    xyz = readdatastream('nodes')
    height = parent.ginput["Geometry"]["height"]
    surf_dist = np.minimum(xyz[:,1], height-xyz[:,1])
    calc_xyz = np.array(composition[0])
    for element in composition[1].keys():
        calc_value = np.array(composition[1][element])
        nodevalues = interp(height-surf_dist, calc_xyz, calc_value) * 100
        nodevalues = np.where(nodevalues < 0, 0, nodevalues)
        adjustdatastream({"Composition_" + element: nodevalues}, datapos="nodes")

def cylinder_interp(parent, composition):
    xyz = readdatastream('nodes')
    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
    calc_xyz = np.array(composition[0])
    for element in composition[1].keys():
        calc_value = np.array(composition[1][element])
        nodevalues = interp(r, calc_xyz, calc_value) * 100
        nodevalues = np.where(nodevalues < 0, 0, nodevalues)
        adjustdatastream({"Composition_" + element: nodevalues}, datapos="nodes")