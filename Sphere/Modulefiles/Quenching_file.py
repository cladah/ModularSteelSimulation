from .ModuleStructure_file import CalcModule


class Quenchingmodule(CalcModule):
    def __init__(self):
        super().__init__("Quenching")

    def run(self):
        from Sphere.Modulefiles.Solvers.ComsolSolver import runComsol

        if self.program == "Comsol":
            print('\nQuenching module')

            if self.program == 'FCSx':
                print('FeniCSx not implemented')
                # rundocker(self)
            elif self.program == 'Comsol':
                print('Using COMSOL for FEM calculation')
                runComsol(self)
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
