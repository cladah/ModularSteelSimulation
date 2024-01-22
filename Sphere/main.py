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
    createdatastreamcache()

    modules = list()
    modules.append(CalcModule("Meshing"))
    modules.append(CalcModule("Carbonitriding"))
    modules.append(CalcModule("TTT"))
    modules.append(CalcModule("TTTmodeling"))
    modules.append(CalcModule("Quenching"))

    for currentmodule in modules:
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=runsinglemodule, args=(currentmodule,))
            tid.start()
            progressmonitor(tid, currentmodule)
        else:
            runsinglemodule(currentmodule)

    removedatastreamcache()

def runsinglemodule(module):
    module.runmodule()

def progressmonitor(tid, module):
    print(module.getprogress())
    if tid.is_alive():
        time.sleep(0.1)
        progressmonitor(tid, module)


if __name__ == "__main__":
    GUI()
