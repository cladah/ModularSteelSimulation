from StructureFile import NewCalcModule


class Carbonitridingmodule(NewCalcModule):
    def __init__(self):
        super().__init__("Carbonitriding")

    def run(self):
        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")
            return

        if self.program == "TC":
            import numpy as np
            from numpy import interp
            from HelpFile import readdatastreamcache, readdatastream, adjustdatastream
            from Solvers.ThermocalcSolver import TCequalibrium, TCcarbonitriding
            print("Carbonitriding module")
            self.updateprogress(0.1)
            if not self.check_runcondition():
                print('Using precalculated carbnitriding simulation')
                for element in self.data["Material"]["Composition"].keys():
                    elementvalues = readdatastreamcache("Composition/" + element)
                    adjustdatastream("Composition/" + element, elementvalues, "nodes")
                print("Carbonitriding module done")
                self.updateprogress(1.0)
                return

            if self.program == "TC":
                print('Running carbon-nitriding module with ThermoCalc')
                activityenv = TCequalibrium("env")
                print("Avtivity of atmosphere calculated")
                self.updateprogress(0.2)
                composition = TCcarbonitriding(activityenv)
                self.updateprogress(0.9)
            else:
                raise KeyError(str(self.program) + ' not implemented for carbonitriding')

            xyz = readdatastream('nodes')
            r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)
            calc_xyz = np.array(composition[0])
            for element in composition[1].keys():
                calc_value = np.array(composition[1][element])
                nodevalues = interp(r, calc_xyz, calc_value) * 100
                adjustdatastream("Composition/" + element, nodevalues, "nodes")
            self.updateprogress(1.0)
            print("Carbonitriding module done")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
