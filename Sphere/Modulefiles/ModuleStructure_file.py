from Framework.Datastream_file import read_input
from Framework.HelpFile import checkruncondition, adjustinputcache


class CalcModule:
    def __init__(self, modulename, runcondition=True, programs=None):
        self.data = read_input()
        self.runcondition = runcondition
        self.module = modulename
        self.progress = 0.0
        self.program = self.data["Programs"][modulename]
        self.__programs = programs

        # Setup methods
        self.check_runcondition()

    def modulename(self):
        return self.module

    def getprogress(self):
        return self.progress

    def updateprogress(self, progressvalue):
        self.progress = progressvalue

    def reset_runcondition(self):
        self.runcondition = True

    def check_runcondition(self):
        self.runcondition = checkruncondition(self.module)
        return self.runcondition

    def __str__(self):
        return self.module

    def possibleprograms(self):
        return self.__programs
