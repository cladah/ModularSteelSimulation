from .ModuleStructure_file_new import CalcModule
import numpy as np
from numpy import interp
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream, getaxisvalues
from Framework.HelpFile import read_input, read_geninput,read_modinput
from .Solvers.ThermocalcSolver import TCequalibrium, TCcarburizing_LPC, TCDiffusionSolver, TC_Cpot, TC_Mob
from .Solvers.FCSxDiff import FCSxDiffSolver

class Diffusionmodule(CalcModule):

    def __init__(self, infile, modulenr=0):
        infile = "Inputs/" + infile + ".json"
        super().__init__("Diffusion", infile, modulenr)

    def run(self):
        Carbtime = sum(self.minput["BoostTime"] + self.minput["DiffTime"]) / 3600
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Diffusion module (General): " + self.inputfile + "\n",
                  "Grainsize is " + str(self.minput["GrainSize"]) + " \u03BCm",
                  "Carburization temperature set to " + str(self.minput["CNtemp"]) + " \N{DEGREE SIGN}C",
                  "Carburization pressure set to " + str(self.minput["CNPress"]) + " kPa",
                  "Carbon potential at surface set to " + str(self.minput["Cpotential"]),
                  str(self.minput["BoostNr"]) + " boost cycles, with boost times of " + str(self.minput["BoostTime"]),
                  "Diffusion times set to " + str(self.minput["DiffTime"]),
                  "Total carburization time is " + str(round(Carbtime, 1)) + " hours",
                  "\n---------------------------------------------------------------------\n\n"]

        for line in outstr:
            self.writeoutput(line)
            print(line)

        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")

            elements = list(self.ginput["Material"]["Composition"].keys())
            #elements.remove("C")
            for element in elements:
                elementvalues = readdatastreamcache("Composition_" + element)
                adjustdatastream({"Composition_" + element: elementvalues}, datapos="nodes")
            print("Carburization module done\n")
            self.updateprogress(1.0)
            return

        if self.program == "TC":
            self.updateprogress(0.1)

            print('Running diffusion module with ThermoCalc')
            activityenv = TCequalibrium(self.ginput, self.minput, "env")
            print("Activity of atmosphere calculated")
            self.updateprogress(0.2)

            composition = dict()
            for element in self.ginput["Material"]["Composition"].keys():
                composition[element] = getaxisvalues("Composition_" + element, -1)


            # Adding geometry to composition?
            composition = TCDiffusionSolver(self.ginput, self.minput, activityenv, composition)
            print(type(composition))
            print("Testing!!!!")

            self.updateprogress(0.9)

        elif self.program == "FCSx":
            self.updateprogress(0.1)

            print('Running diffusion module with ThermoCalc and FeniCSx')
            mobC = TC_Mob(self.ginput, self.minput)

            print("Activity of atmosphere calculated")
            self.updateprogress(0.2)

            composition = dict()
            for element in self.ginput["Material"]["Composition"].keys():
                composition[element] = readdatastream("Composition_" + element, -1)

            # Adding geometry to composition?
            composition = FCSxDiffSolver(self.ginput, self.minput, mobC, composition)
            self.updateprogress(0.9)

            for element in composition.keys():
                adjustdatastream({"Composition_" + element: composition[element]}, datapos="nodes")
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
        if self.ginput["Geometry"]["Type"] == "4PointBend":
            fourpointbend_interp(self, composition)
        elif self.ginput["Geometry"]["Type"] == "2Daxisym":
            cylinder_interp(self, composition)
        else:
            raise KeyError("Type of geometry not supported for compositional "
                           "interpolation in diffusion modell")

        self.updateprogress(1.0)
        print("Carburization module done\n")
def fourpointbend_interp(parent, composition):
    print("interpolation testing!")
    xyz = readdatastream('nodes')
    height = parent.ginput["Geometry"]["height"]
    surf_dist = np.minimum(xyz[:,1], height-xyz[:,1])
    print(np.min(surf_dist))
    print(np.max(surf_dist))
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