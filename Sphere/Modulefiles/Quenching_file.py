from .ModuleStructure_file import CalcModule


class Quenchingmodule(CalcModule):
    def __init__(self):
        super().__init__("Quenching")

    def run(self):
        from Sphere.Modulefiles.Solvers.ComsolSolver import runComsol

        if self.program == "Comsol":
            runComsol(self)
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
