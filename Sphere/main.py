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
        print(" ")
        if programstate.get() == 0 or runall == 1:
            createMesh()
            addgui(gui, "Mesh")
        elif programstate.get() == 1 or runall == 1:
            runcarbonitridingmodule()
            addgui(gui, "Carbonitriding")
        elif programstate.get() == 2 or runall == 1:
            runTTTmodule()
            addgui(gui, "TTT")
        elif programstate.get() == 3 or runall == 1:
            runTTTmodelmodule()
            addgui(gui, "TTTmodel")
        elif programstate.get() == 4 or runall == 1:
            runquenchingmodule()
            addgui(gui, "Quenching")
        programstate.set(programstate.get()+1)
    if not os.path.exists('Cachefiles/InputCache.json'):
        createinputcache()
    gui = tk.Tk()
    gui.geometry("500x600")
    gui.title("Quenching of steel")
    gui = runguimodule(gui)
    button = tk.Button(master=gui, text="Continue", command=runmodule, padx=20)
    button.pack(side=tk.TOP, anchor=tk.NW)
    programstate = tk.IntVar(gui, 0)
    runall = tk.IntVar(gui, 0)

    gui.mainloop()

def modelling():
    createMesh()
    runcarbonitridingmodule()
    runTTTmodule()
    runTTTmodelmodule()
    runquenchingmodule()

if __name__ == "__main__":
    start()
