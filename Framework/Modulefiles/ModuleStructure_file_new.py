from Framework.HelpFile import checkruncondition, adjustinputcache, read_geninput, read_modinput, print_output
import pathlib

class CalcModule:
    def __init__(self, modulename, inputfile, modulenr=0, runcondition=True):
        self.ginput = read_geninput()
        self.minput = read_modinput(inputfile)
        self.runcondition = runcondition
        self.module = modulename
        self.progress = 0.0
        self.program = self.ginput["Programs"][modulename]
        self.inputfile = inputfile
        #file = pathlib.Path(inputfile).name.removesuffix((".json"))
        #self.modulenr = self.ginput["Inputs"].index(file)
        self.modulenr = modulenr
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

    def __repr__(self):
        return "<" + str(self.module) + " module object>"
    def __str__(self):
        return "<" + str(self.module) + " module object>"

    def writeoutput(self, outputstr):
        print_output(outputstr)
