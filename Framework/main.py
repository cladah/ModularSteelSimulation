import time
import tkinter

import meshio
import numpy as np
from matplotlib import pyplot as plt

from CGUImodule import MainApp
from Result_GUI import Result_MainApp, Compare_MainApp
from Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input, createinputcache, change_input, reset_input, analyseTTTdatabase, get_plotlbls, read_geninput, reset_output
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule, Carbonizationmodule, Diffusionmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule
from Framework.Modulefiles.Docker_file import rundocker, rundocker_1D
from Modulefiles.Testmod_file import TestModule
from Datastream_file import getaxisvalues, readdatastream
from ResultReading import read_results_history, read_results, getnames_results, read_results_all, read_results_axis
from Postprocessing.dataextraction import DatastreamPlotting, ResultPlotting, export_data

def progressmonitor(tid, module):
    """
    :param tid:
    :param module:
    """
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)


def GUI():
    """
    Using Tkinter to run the GUI
    """
    print("Opening GUI window...")
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()

def Result_GUI_comp(filenames, dataname):
    """
            Using Tkinter to run the GUI
        """
    print("Opening result window...")
    ctk.set_appearance_mode("dark")
    app = Compare_MainApp(filenames, dataname)
    app.mainloop()

def Result_GUI_show(filename):
    """
        Using Tkinter to run the GUI
    """
    if filename == "":
        root = tkinter.Tk()
        root.withdraw()
        filename = tkinter.filedialog.askopenfilename(filetypes=(("xdmf files", "*.xdmf"),))
    print("Opening result window...")
    ctk.set_appearance_mode("dark")
    app = Result_MainApp(filename)
    app.mainloop()


def looping():
    """
        Running multiple simulations with specified modules

        Appending modules are done by hand.

        :return: xdmf file with simulation results
        """

    setupSimulation()

    differentin = [["Thermo", "CNtemp", 1173.15], ["Thermo", "CNtemp", 973.15], ["Thermo", "CNtime", 172800]]
    saveloc = ["October2024_900C.xdmf", "October2024_700C.xdmf", "October2024_CN2Days.xdmf"]

    if len(differentin) != len(saveloc):
        raise KeyError("Wrong looping input")

    i = 0
    for a in differentin:
        print("Loop nr " + str(i+1))
        reset_input("Cachefiles/Input_ref.json")
        change_input(*a)
        data = read_input()
        createdatastreamcache(data["Datastream"]["Cachedirect"])

        # resetdatastream()
        createinputcache()

        modules = list()
        modules.append(Meshingmodule())
        modules.append(Carbonizationmodule())
        # modules.append(TTTdiagrammodule())
        # modules.append(Transformationmodelmodule())
        #modules.append(Quenchingmodule())

        for currentmodule in modules:
            currentmodule.run()

        removedatastreamcache()
        savedatastream(saveloc[i])
        i = i + 1


def modelling():
    """
    Running a single simulation with specified modules

    Appending modules are done by hand.

    :return: xdmf file with simulation results
    """

    setupSimulation()

    data = read_input()


    # resetdatastream()
    createinputcache()

    modules = list()
    modules.append(Meshingmodule())
    modules.append(Carbonizationmodule())
    modules.append(TTTdiagrammodule())
    modules.append(Transformationmodelmodule())
    modules.append(Quenchingmodule())

    for currentmodule in modules:
        currentmodule.run()
        savedatastream(data["Datastream"]["Savedirect"])

    removedatastreamcache()
    savedatastream(data["Datastream"]["Savedirect"])
    Result_GUI_show("Datastream.xdmf")

def DockerTest():
    print("Running Docker testing env")
    #rundocker()
    rundocker_1D()
    pass
def vtxfile_test():
    from adios2 import FileReader

    with FileReader("FeniCSx/beam_stress.bp") as s:
        # inspect variables
        vars = s.available_variables()
        for name, info in vars.items():
            print("variable_name: " + name, end=" ")
            for key, value in info.items():
                print("\t" + key + ": " + value, end=" ")
            print()
        print()
    pass

def FCSxfile_test():
    with meshio.xdmf.TimeSeriesReader("FeniCSx/beam_stress.xdmf") as reader:
        #points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            print(k)
            print(reader.domain)
            #t, point_data, cell_data = reader.read_data(k)
            #datanames = point_data.keys()
            #print(datanames)

def setupSimulation():
    """
    Simulation setup from information in input.json

    :return:
    modules - Modules in order of exicution (list)
    """

    print("Setting up simulation...")
    reset_output()
    ginput = read_geninput()
    createdatastreamcache(ginput["Datastream"]["Cachedirect"])
    createinputcache()
    modules = list()
    for i in range(len(ginput["Modules"])):
        infile = ginput["InputDirectory"] + "/" + ginput["Inputs"][i]
        if ginput["Modules"][i] == "Meshing":
            modules.append(Meshingmodule(infile))
        elif ginput["Modules"][i] == "Test":
            modules.append(TestModule(infile))
        elif ginput["Modules"][i] == "Diffusion":
            modules.append(Diffusionmodule(infile))
        elif ginput["Modules"][i] == "TTTdiagram":
            modules.append(TTTdiagrammodule(infile))
        elif ginput["Modules"][i] == "TransformMod":
            modules.append(Transformationmodelmodule(infile))
        elif ginput["Modules"][i] == "Quenching":
            modules.append(Quenchingmodule(infile))
        else:
            raise KeyError("Module input in iMain not supported")
    print("Simulation structure setup.")
    return modules

if __name__ == "__main__":
    """
    ginput - General input (dict)
    """
    run = input("What do you want to run? 1 - Run 2 - GUI 3 - Result\n")
    print(run)
    print(type(run))
    if int(run) == 1:
        print("Running normal input")
        ginput = read_geninput()
        modules = setupSimulation()
        for module in modules:
            module.run()
            savedatastream(ginput["Datastream"]["Savedirect"])
        Result_GUI_show("Resultfiles/" + ginput["Datastream"]["Savedirect"])
    if run == 2:
        GUI()
        #DockerTest()
    elif run == 3:
        Result_GUI_show("")
        #export_data("Resultfiles/159A_Carb3.xdmf", ["Composition/C", "Martensite"], -1)
        #export_data("Resultfiles/159A_Carb3.xdmf", "All", -1)
    if False:
        files = ["Resultfiles/October2024_LPC_4h_2.xdmf", "Resultfiles/October2024_LPC_2h.xdmf"]
        dataname = "Composition/C"
        ResultPlotting(files, dataname)

