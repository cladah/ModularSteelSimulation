import threading


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
        self.progress = 0.0
    def modulename(self):
        return self.module
    def getprogress(self):
        return self.progress
    def updateprogress(self, progressvalue):
        self.progress = progressvalue
    def reset(self):
        self.run = True
        pass

    def runmodule(self):
        self.updateprogress(0.0)
        if self.module == "Test":
            print("Testing module")
            self.testingmodule()
        elif self.module == "Meshing":
            from Meshmodule import createMesh
            createMesh(self)
        elif self.module == "Carbonitriding":
            from Carbonitridingmodule import runcarbonitridingmodule
            runcarbonitridingmodule(self)
        elif self.module == "TTT":
            from TTTmodule import runTTTmodule
            runTTTmodule(self)
        elif self.module == "TTTmodeling":
            from TTTmodule import runTTTmodelmodule
            runTTTmodelmodule(self)
        elif self.module == "Quenching":
            from Quenchingmodule import runquenchingmodule
            runquenchingmodule(self)
        else:
            print(self.module + " not implemented.")
        self.updateprogress(1.0)
    def testingmodule(self):
        import time
        for i in range(1, 3):
            time.sleep(5)
            print(str(i * 5) + "sec")