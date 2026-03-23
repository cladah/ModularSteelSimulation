from Framework.HelpFile import checkruncondition, adjustinputcache, read_geninput, read_modinput, print_output
import pathlib
from typing import Any

class CalcModule_old:
    def __init__(self, module_name: str, input_file: str, module_nr: int = 0, run_condition: bool = True):
        # Using underscores for attribute naming (PEP 8 standard)
        self.ginput = read_geninput()
        self.minput = read_modinput(input_file)
        self.module = module_name
        self.input_file = input_file
        self.module_nr = module_nr

        # Initialize internal state
        self.program = self.ginput.get("Programs", {}).get(module_name)

        # Initial run condition check
        self.run_condition = run_condition
        self.check_run_condition()

    @property
    def progress(self) -> float:
        """Getter for progress."""
        return self._progress

    @progress.setter
    def progress(self, value: float):
        """Setter for progress with basic validation."""
        if 0.0 <= value <= 100.0:
            self._progress = value
        else:
            self._progress = value  # Or raise a ValueError if strict

    def reset_run_condition(self):
        """Resets the run condition to True."""
        self.run_condition = True

    def check_run_condition(self) -> bool:
        """Updates and returns the current run condition."""
        self.run_condition = checkruncondition(self.module)
        return self.run_condition

    def write_output(self, output_str: str):
        print_output(output_str)




class CalcModule:
    def __init__(self, modulename, inputfile, modulenr=0, runcondition=True):
        self.ginput = read_geninput()
        self.minput = read_modinput(inputfile)
        self.runcondition = runcondition
        self.module = modulename
        self._progress = 0.0
        self.program = self.ginput["Programs"][modulename]
        self.inputfile = inputfile
        self.modulenr = modulenr
        self.check_runcondition()
        #self.program = self.ginput.get("Programs", {}).get(modulename)

    @property
    def progress(self) -> float:
        return self._progress

    @progress.setter
    def progress(self, value: float):
        if 0.0 <= value <= 100.0:
            self._progress = value
        else:
            self._progress = value

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

    def __repr__(self) -> str:
        return f"<CalcModule(name='{self.module}', nr={self.module_nr})>"

    def __str__(self) -> str:
        return f"Module: {self.module} (Progress: {self.progress}%)"

    def writeoutput(self, outputstr):
        print_output(outputstr)
