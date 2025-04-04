from .ModuleStructure_file_new import CalcModule
from .Solvers.ThermocalcSolver import getTTTcompositions
from Framework.HelpFile import getTTTdata, addTTTdata
from .Solvers.ThermocalcSolver import calculatePearlite, calculateBainite, calculateFerrite, calculateMartensite
import numpy as np


class TTTdiagrammodule(CalcModule):
    def __init__(self, infile):
        infile = "Inputs/" + infile + ".json"
        super().__init__("TTTdiagram", infile)

    def run(self):
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Transformation diagram module: " + self.inputfile + "\n",
                  "Grainsize is " + str(self.minput["GrainSize"]) + " \u03BCm",
                  "Upper temperature set to: " + str(1000) + " \N{DEGREE SIGN}C",
                  "Lower temperature set to: " + str(0) + " \N{DEGREE SIGN}C",
                  "Number of steps: " + str(40),
                  "\n---------------------------------------------------------------------\n\n"]

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
                runTTTcalc(tmpcomp)
                self.updateprogress(i / compnr)
                i = i + 1
            print("TTT diagram module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

def runTTTcalc(composition):
    # Take away values that are 0 EXCEPT N AND C!
    tmpcomp = composition.copy()
    keys = tmpcomp.keys()
    for key in keys:
        if key in ["C", "N"]:
            continue
        if tmpcomp[key] == 0.0:
            del composition[key]
    #
    if bool(getTTTdata(composition, "TTTdata")):
        print("TTTdata exists in database for " + str(composition))
        return
    print("Running TTT calculation for " + str(composition))
    # If N is changed check
    tmpcomp = composition.copy()
    tmpcomp["N"] = 0.
    if bool(getTTTdata(tmpcomp, "TTTdata")):
        TTTdata = getTTTdata(tmpcomp, "TTTdata")
        start, half, finish = calculateMartensite(composition)
        start = [[start, start], [0.1, 1E12]]
        half = [[half, half], [0.1, 1E12]]
        finish = [[finish, finish], [0.1, 1E12]]
        phase = dict()
        phase["start"] = start
        phase["half"] = half
        phase["finish"] = finish
        #print(TTTdata)
        #print(phase)
        TTTdata["Martensite"] = phase

        addTTTdata(composition, TTTdata, "TTTdata")
        return

    Tsteps = np.linspace(0, 1000, 41) + 273.15
    phases = ["Ferrite", "Bainite", "Pearlite", "Martensite"]
    TTTdata = dict()
    for ph in phases:
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
        phase = dict()
        phase["start"] = start
        phase["half"] = half
        phase["finish"] = finish
        TTTdata[ph] = phase
    addTTTdata(composition, TTTdata, "TTTdata")