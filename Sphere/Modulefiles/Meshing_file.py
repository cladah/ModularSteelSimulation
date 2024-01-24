from StructureFile import NewCalcModule


class Meshingmodule(NewCalcModule):
    def __init__(self):
        super().__init__("Meshing")

    def run(self):
        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")
            return

        from Solvers.MeshSolvers import gmshsolver, pygmshsolver

        if self.program == "Gmsh":
            gmshsolver(self)
        elif self.program == "PyGmsh":
            pygmshsolver(self)
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
