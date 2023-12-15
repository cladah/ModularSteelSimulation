import os

import numpy as np
import h5py
import matplotlib.pyplot as plt

from HelpFile import *
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitridingmodule
from TTTmodule import runTTTmodule, runTTTmodelmodule
from GUImodule import runguimodule, addgui
from Postmodule import *
from Quenchingmodule import runquenchingmodule
from Testfile import *
import tkinter as tk



def start():
    def runmodule():
        print(programstate.get())
        if programstate.get() == 0:
            createMesh()
            addgui(gui, "Mesh")
        elif programstate.get() == 1:
            runcarbonitridingmodule()
            print(getaxisvalues("Composition/C"))
            input("")
        elif programstate.get() == 2:
            runTTTmodule()
        elif programstate.get() == 3:
            runTTTmodelmodule()
        elif programstate.get() == 4:
            modeldatatocsv()
        elif programstate.get() == 5:
            runquenchingmodule()
        programstate.set(programstate.get()+1)
    if not os.path.exists('Cachefiles/InputCache.json'):
        createinputcache()
    gui = runguimodule()
    button = tk.Button(master=gui, text="Run", command=runmodule, padx=20)
    button.pack(side=tk.TOP, anchor=tk.SW)
    programstate = tk.IntVar(gui, 0)
    gui.mainloop()

def modelling():
    createMesh()
    runcarbonitridingmodule()
    runTTTmodule()
    runTTTmodelmodule()
    modeldatatocsv()
    runquenchingmodule()

if __name__ == "__main__":
    modelling()
