from .ModuleStructure_file_new import CalcModule
from .Solvers.ThermocalcSolver import getTTTcompositions
from Framework.HelpFile import getTTTdb, addTTTdb, checkTTTdb
from .Solvers.ThermocalcSolver import calculatePearlite, calculateBainite, calculateFerrite, calculateMartensite
import numpy as np


class TTTdiagrammodule(CalcModule):
    def __init__(self, infile, modulenr):
        infile = "Inputs/" + infile + ".json"
        super().__init__("TTTdiagram", infile, modulenr)

    def run(self):
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Transformation diagram module: " + self.inputfile + "\n",
                  "Grainsize is " + str(self.minput["GrainSize"]) + " \u03BCm",
                  "Upper temperature set to: " + str(1000) + " \N{DEGREE SIGN}C",
                  "Lower temperature set to: " + str(0) + " \N{DEGREE SIGN}C",
                  "Number of steps: " + str(40)]

        for ph in self.minput["Phases"]:
            growth_mode = self.minput["GrowthMode"]
            outstr.append(f"{ph} modeled with {growth_mode[ph]}")

        outstr.append("\n---------------------------------------------------------------------\n")

        for line in outstr:
            self.writeoutput(line)
            print(line)

        if not self.runcondition:
            print("Using precalculated " + str(self.module) + " simulation")
            print("TTT diagram module done\n")
            return

        if self.program == "TC":
            print("\nTTT module")
            TTTcompositions = getTTTcompositions()
            compnr = len(TTTcompositions)
            i = 1
            print("Number of TTT calculations are " + str(compnr))
            for tmpcomp in TTTcompositions:
                runTTTcalc(tmpcomp, self.minput)
                self.updateprogress(i / compnr)
                i = i + 1
            print("TTT diagram module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

def runTTTcalc(composition, minput):
    # Take away values that are 0 EXCEPT N AND C!
    tmpcomp = composition.copy()
    keys = tmpcomp.keys()
    for key in keys:
        if key in ["C", "N"]:
            continue
        if tmpcomp[key] == 0.0:
            del composition[key]

    model_input = {k: minput[k] for k in ["GrainSize"]}

    # If N is changed check
    tmpcomp = composition.copy()
    tmpcomp["N"] = 0.

    Tsteps = np.linspace(0, 1000, 41) + 273.15
    phases = list(minput["Phases"])

    print("Running TTT calculation for " + str(composition))
    for ph in phases:
        growth_mode = minput["GrowthMode"][ph]
        if checkTTTdb(composition, model_input, f"{ph}_{growth_mode}"):
            continue
        if ph == "Ferrite":
            start, half, finish = calculateFerrite(Tsteps, composition)
            start = [Tsteps, start]
            half = [Tsteps, half]
            finish = [Tsteps, finish]
        elif ph == "Pearlite":
            start, half, finish = calculatePearlite(Tsteps, composition)
            start = [Tsteps, start]
            half = [Tsteps, half]
            finish = [Tsteps, finish]
        elif ph == "Bainite":
            start, half, finish = calculateBainite(Tsteps, composition)
            start = [Tsteps, start]
            half = [Tsteps, half]
            finish = [Tsteps, finish]
        elif ph == "Martensite":
            start, half, finish = calculateMartensite(composition)
            start = [[start, start], [0.1, 1E12]]
            half = [[half, half], [0.1, 1E12]]
            finish = [[finish, finish], [0.1, 1E12]]
        else:
            raise KeyError("Phase error in TTT module, check input phases")
        data = dict()
        data["start"] = start
        data["half"] = half
        data["finish"] = finish

        addTTTdb(data, composition, model_input, f"{ph}_{growth_mode}")
    print("\n")