from .ModuleStructure_file import CalcModule


class Transformationmodelmodule(CalcModule):
    def __init__(self):
        super().__init__("Quenching")

    def run(self):
        if not self.runcondition:
            print("Using precalculated " + str(self.module) + " simulation")
            return

        if self.program == "Python":
            print()
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
