from .ModuleStructure_file_new import CalcModule
import numpy as np
from numpy import interp
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream
from Framework.HelpFile import read_input
from .Solvers.ThermocalcSolver import TCequalibrium, TCcarbonitriding, TCcarburizing, TCcarburizing_LPC

class Carbonizationmodule(CalcModule):
    def __init__(self):
        super().__init__("Carburization", "Cachefiles/iDiff.json")

    def run(self):
        Carbtime = sum(self.minput["BoostTime"] + self.minput["DiffTime"])/3600
        outstr = ["Diffusion module",
                  "Grainsize is " + str(self.minput["GrainSize"]) + " \u03BCm",
                  "Carburization temperature set to " + str(self.minput["CNtemp"]) + " \N{DEGREE SIGN}C",
                  "Carburization pressure set to " + str(self.minput["CNPress"]) + " kPa",
                  "Carbon potential set to " + str(self.minput["Cpotential"]),
                  str(self.minput["BoostNr"]) + " boost cycles, with boost times of " + str(self.minput["BoostTime"]),
                  "Diffusion times set to " + str(self.minput["DiffTime"]),
                  "Total carburization time is " + str(round(Carbtime, 1)) + " hours",
                  "---------------------------------------------------------------------\n\n"]

        for line in outstr:
            self.writeoutput(line)
            print(line)


        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")
            for element in self.data["Material"]["Composition"].keys():
                elementvalues = readdatastreamcache("Composition/" + element)
                adjustdatastream({"Composition/" + element: elementvalues}, "nodes")
            print("Carburization module done\n")
            self.updateprogress(1.0)
            return

        if self.program == "TC":
            print("Carburization module")
            data = read_input()
            self.updateprogress(0.1)

            print('Running carburization module with ThermoCalc')
            activityenv = TCequalibrium("env")
            print("Avtivity of atmosphere calculated")
            self.updateprogress(0.2)

            """
                Choosing which carburizing model to run.
                
                P<0.01atm -> Low pressure carurizing
                P>0.01atm -> Atmosphere pressure carurizing
            """
            if data["Thermo"]["CNPress"] > 10000:
                composition = TCcarburizing(activityenv)
            else:
                boosts = data["Thermo"]["BoostNr"]
                composition = TCcarburizing_LPC(activityenv, boosts, data["Thermo"]["BoostTime"], data["Thermo"]["DiffTime"])
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

class Carbonitridingmodule(CalcModule):
    def __init__(self):
        super().__init__("Carbonitriding")

    def run(self):
        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")
            for element in self.data["Material"]["Composition"].keys():
                elementvalues = readdatastreamcache("Composition/" + element)
                adjustdatastream({"Composition/" + element: elementvalues}, "nodes")
            print("Carbonitriding module done\n")
            self.updateprogress(1.0)
            return

        if self.program == "TC":
            print("Carbonitriding module")
            self.updateprogress(0.1)

            print('Running carbon-nitriding module with ThermoCalc')
            activityenv = TCequalibrium("env")
            print("Activity of atmosphere")
            print(activityenv)
            print("Activity of atmosphere calculated")
            self.updateprogress(0.2)

            """
                Choosing which carbonitriding model to run.

                P<0.01atm -> Low pressure carbonitriding
                P>0.01atm -> Atmosphere pressure carbonitriding
            """

            composition = TCcarbonitriding(activityenv)
            self.updateprogress(0.9)

            xyz = readdatastream('nodes')
            r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
            calc_xyz = np.array(composition[0])
            for element in composition[1].keys():
                calc_value = np.array(composition[1][element])
                nodevalues = interp(r, calc_xyz, calc_value) * 100
                nodevalues = np.where(nodevalues < 0, 0, nodevalues)
                adjustdatastream({"Composition/" + element: nodevalues}, "nodes")
            self.updateprogress(1.0)
            print("Carbonitriding module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
