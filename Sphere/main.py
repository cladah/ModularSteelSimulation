import os
import time

import matplotlib.pyplot as plt
from CGUImodule import MainApp

from Testfile import *
import customtkinter as ctk
from StructureFile import CalcModule


def GUI():
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()
def modelling():
    import threading
    createdatastreamcache("Resultfiles/TestDatastream.xdmf")

    modules = list()
    modules.append(CalcModule("Meshing"))
    modules.append(CalcModule("Carbonitriding"))
    modules.append(CalcModule("TTT"))
    modules.append(CalcModule("TTTmodeling"))
    modules.append(CalcModule("Quenching"))

    for currentmodule in modules:
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=run_single_module, args=(currentmodule,))
            tid.start()
            progressmonitor(tid, currentmodule)
        else:
            run_single_module(currentmodule)
    input("")
    removedatastreamcache()
    savedatastream("Resultfiles/TestDatastream.xdmf")

def run_single_module(module):
    module.runmodule()

def progressmonitor(tid, module):
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)
if __name__ == "__main__":
    GUI()
