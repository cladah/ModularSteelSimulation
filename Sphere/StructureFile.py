
class CalcModule:
    def __init__(self, moduletype, run=True):
        from HelpFile import read_input
        self.data = read_input()
        self.module = moduletype
        try:
            self.program = self.data["Programs"][moduletype]
        except KeyError:
            print(list(self.data["Programs"].keys()))
            raise KeyError("Moduletype only takes arguments " + str(list(self.data["Programs"].keys())))
        self.run = run

    def reset(self):
        self.run = True
        pass

    def runmodule(self):
        pass
        if self.module == "Test":
            import time
            print("Testing module")
            for i in range(1, 3):
                time.sleep(5)
                print(str(i*5)+"sec")
        elif self.module == "Meshing":
            from Meshmodule import createMesh
            createMesh()
        elif self.module == "Carbonitriding":
            from Carbonitridingmodule import runcarbonitridingmodule
            runcarbonitridingmodule()
        elif self.module == "TTT":
            from TTTmodule import runTTTmodule
            runTTTmodule()
        elif self.module == "TTTmodeling":
            from TTTmodule import runTTTmodelmodule
            runTTTmodelmodule()
        elif self.module == "Quenching":
            from Quenchingmodule import runquenchingmodule
            runquenchingmodule()
        else:
            print(self.module + " not implemented.")