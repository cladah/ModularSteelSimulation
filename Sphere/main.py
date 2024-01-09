import os

import numpy as np
import h5py
import matplotlib.pyplot as plt

from HelpFile import *
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitridingmodule
from TTTmodule import runTTTmodule, runTTTmodelmodule
from GUImodule import runguimodule, addgui, MainApp
from Postmodule import *
from Quenchingmodule import runquenchingmodule
from Testfile import *
import tkinter as tk
from tkinter import ttk



def start():
    data = read_input()
    def runmodule():
        print(" ")
        if programstate.get() == 0 or runall == 1:
            createdatastreamcache()
            createMesh()
            addgui(tabs, "Mesh")
            tabs.select(1)
        elif programstate.get() == 1 or runall == 1:
            runcarbonitridingmodule()
            addgui(tabs, "Carbonitriding")
            tabs.select(2)
        elif programstate.get() == 2 or runall == 1:
            runTTTmodule()
            addgui(tabs, "TTT")
            tabs.select(3)
        elif programstate.get() == 3 or runall == 1:
            runTTTmodelmodule()
            addgui(tabs, "TTTmodel")
            tabs.select(4)
        elif programstate.get() == 4 or runall == 1:
            runquenchingmodule()
            removedatastreamcache()
            addgui(tabs, "Quenching")
            tabs.select(5)
            button.destroy()
        elif programstate.get() > 4 or runall == 1:
            print("All simulations are done")
        programstate.set(programstate.get()+1)
    if not os.path.exists('Cachefiles/InputCache.json'):
        createinputcache()
    gui = tk.Tk()
    gui.geometry("1200x1000")
    gui.title("Quenching of steel")
    header = tk.Frame(gui, height=10)
    button = tk.Button(master=header, text="Continue", command=runmodule, padx=20, pady=10)
    button.pack(side=tk.TOP, anchor=tk.NW, padx=20,pady=20)
    composition = tk.Label(header, text="Composition of steel is:", pady=20)
    composition.pack()
    compstr = " ".join(
        [i + "=" + str(data["Material"]["Composition"][i]) for i in data["Material"]["Composition"].keys()])
    composition = tk.Label(header, text=compstr)
    composition.pack()


    header.pack(side="top", fill="both", expand=False)

    # Add
    tabs = ttk.Notebook(gui)
    tabs.pack(side="top", fill="both", expand=True)

    runguimodule(tabs)

    #container = tk.Frame(gui)
    #container.pack(side="top", fill="both", expand=True)
    #container.grid_rowconfigure(0, weight=1)
    #container.grid_columnconfigure(0, weight=1)

    programstate = tk.IntVar(gui, 0)
    runall = tk.IntVar(gui, 0)

    gui.mainloop()
def GUI():
    app = MainApp()
    app.mainloop()
def modelling():
    createdatastreamcache()
    createMesh()
    runcarbonitridingmodule()
    runTTTmodule()
    runTTTmodelmodule()
    runquenchingmodule()
    removedatastreamcache()

if __name__ == "__main__":
    start()
