from .ModuleStructure_file import CalcModule
from .Solvers.ThermocalcSolver import getTTTcompositions
from Framework.HelpFile import getTTTdata, addTTTdata
from .Solvers.ThermocalcSolver import calculatePearlite, calculateBainite, calculateFerrite, calculateMartensite
import numpy as np


class TTTdiagrammodule(CalcModule):
    def __init__(self):
        super().__init__("TTT")

    def run(self):
        if not self.runcondition:
            print("Using precalculated " + str(self.module) + " simulation")
            print("TTT diagram module done\n")
            return

        if self.program == "TC":
            print("\nTTT module")
            TTTcompositions = getTTTcompositions()
            compnr = len(TTTcompositions)
            i = 1
            for tmpcomp in TTTcompositions:
                print(tmpcomp)
                runTTTcalc(tmpcomp)
                self.updateprogress(i / compnr)
                i = i + 1
            print("TTT diagram module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

def runTTTcalc(composition):
    # Take away values that are 0 EXCEPT N AND C!
    keys = composition.keys()
    for key in keys:
        if key in ["C", "N"]:
            continue
        if composition[key] == 0.:
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