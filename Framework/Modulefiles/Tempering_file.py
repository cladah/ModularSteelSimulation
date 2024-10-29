from Framework.Oldfiles.StructureFile import NewCalcModule


class Temperingmodule(NewCalcModule):
    def __init__(self):
        super().__init__("Tempering")

    def run(self):

        if self.program == "":
            pass
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
