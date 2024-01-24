from StructureFile import NewCalcModule


class Transformationmodelmodule(NewCalcModule):
    def __init__(self):
        super().__init__("Quenching")

    def run(self):
        if not self.runcondition:
            print("Using precalculated " + str(self.module) + " simulation")
            return

        from ..Solvers.ComsolSolver import runComsol

        if self.program == "Python":
            print()
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
