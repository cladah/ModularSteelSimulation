from .ModuleStructure_file_new import CalcModule
import numpy as np
from numpy import interp
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream, getaxisvalues
from Framework.HelpFile import read_input
from .Solvers.ThermocalcSolver import TCequalibrium, TCcarburizing_LPC, TCDiffusionSolver

class Diffusionmodule(CalcModule):

    def __init__(self, infile):
        infile = "Inputs/" + infile + ".json"
        super().__init__("Diffusion", infile)

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
            for element in self.ginput["Material"]["Composition"].keys():
                elementvalues = readdatastreamcache("Composition/" + element)
                adjustdatastream({"Composition/" + element: elementvalues}, "nodes")
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
                composition[element] = getaxisvalues("Composition/" + element, -1)
            # Adding geometry to composition?
            composition = TCDiffusionSolver(self.ginput, self.minput, activityenv, composition)
            self.updateprogress(0.9)

            """
                Interpolating the compositional values to nodal points 
            """

            xyz = readdatastream('nodes')
            r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
            calc_xyz = np.array(composition[0])
            for element in composition[1].keys():
                calc_value = np.array(composition[1][element])
                nodevalues = interp(r, calc_xyz, calc_value) * 100
                nodevalues = np.where(nodevalues < 0, 0, nodevalues)
                adjustdatastream({"Composition/" + element: nodevalues}, "nodes")
            self.updateprogress(1.0)
            print("Carburization module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

