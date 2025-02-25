from .ModuleStructure_file_new import CalcModule

class TestModule(CalcModule):
    def __init__(self):
        super().__init__("Test", "Cachefiles/iDiff_1.json")

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
                  "\n---------------------------------------------------------------------\n\n"]

        for line in outstr:
            self.writeoutput(line)
            print(line)

        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")
            self.updateprogress(1.0)
            return

