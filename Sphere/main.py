import time
import meshio
import numpy as np
from matplotlib import pyplot as plt

from CGUImodule import MainApp
from Result_GUI import Result_MainApp
from Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input, setupSimulation, createinputcache, change_input, reset_input
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule, Carbonizationmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule
from Sphere.Modulefiles.Docker_file import rundocker

from Datastream_file import getaxisvalues, readdatastream
from ResultReading import read_results_history, read_results, getnames_results, read_results_all, read_results_axis


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


def Result_GUI_show(filename):
    """
        Using Tkinter to run the GUI
    """
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
    """
    differentin = [
        ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.3, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.4, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.5, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.4, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.3, "Mo": 0.3, "Si": 0.2}]]  # Double htc

    saveloc = ["Ref.xdmf", "C03.xdmf", "C04.xdmf", "C05.xdmf", "Cr14.xdmf", "Ni13.xdmf"]
    """

    differentin = [["CNtime", 86400], ["CNtime", 86400*2]]


    differentin = [["Thermo", "CNtemp", 1173.15], ["Thermo", "CNtemp", 973.15]]
    saveloc = ["October2024_900C.xdmf", "October2024_700C.xdmf"]

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
        modules.append(TTTdiagrammodule())
        modules.append(Transformationmodelmodule())
        modules.append(Quenchingmodule())

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

    createdatastreamcache(data["Datastream"]["Cachedirect"])
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

    removedatastreamcache()
    savedatastream(data["Datastream"]["Savedirect"])


def DockerTest():
    print("Running Docker testing env")
    mod = Meshingmodule()
    mod.run()
    rundocker()
    pass


def xmdftesting():

    t_data = 0

    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()
        pd_list, cd_list, t_list = list(), list(), list()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            t_list.append(t)
            pd_list.append(point_data)
            cd_list.append(cell_data)
    print(len(pd_list[0]["Composition/C"]))
    print(np.array(np.ones((len(pd_list[0]["Composition/C"]), 5))))

    data = {"ones":np.array(np.ones((len(pd_list[0]["Composition/C"]), 5)))}
    data = {"ones":np.array([[1,1,1,1,1,1] for i in range(len(pd_list[0]["Composition/C"]))])}

    if t_data not in t_list:
        t_list.append(t_data)
        pd_list.append(data)
        cd_list.append({})
    else:
        indx = t_list.index(t_data)
        for key in data.keys():
            pd_list[indx][key] = data[key]

    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(points, cells)
        for i in range(len(t_list)):
            writer.write_data(t=t_list[i], point_data=pd_list[i], cell_data=cd_list[i])


def ResultfileTest():
    import matplotlib.pyplot as plt
    import csv
    import pandas as pd
    def smoothen(x, winsize=5):
        return np.array(pd.Series(x).rolling(winsize).mean())[winsize - 1:]

    with open("Resultfiles/Testing.txt", 'r') as file:
        #reader = csv.reader(file, delimiter="\t")
        r = []
        y = []
        for row in file.readlines():
            row = row.split()
            r.append(np.sqrt(float(row[0])**2+float(row[1])**2))
            y.append(float(row[2]))
        r, y = zip(*sorted(zip(r, y)))
        plt.plot(smoothen(r,50), smoothen(y,50))
        plt.show()


def ResultPlotting(filenames, dataname, point=[0., 0.], tid = -1):
    import matplotlib.pyplot as plt
    data = list()
    for filename in filenames:
        names = getnames_results(filename)
        print(filename)
        print(names)
        tmpdata = read_results_axis(filename, dataname, tid)
        xyz = read_results_axis(filename, "nodes")
        plt.plot(xyz[:,0], tmpdata)
    plt.show()
    print("Done")
    return
    for y in Y:
        data_t, data = read_results_history(filename, dataname, y)
        leg.append(str(round(y[0], 4)))
        plt.plot(data_t, np.array(data)[:,0])
        #plt.plot(data_t, np.array(data)[:,0])
    plt.legend(leg)
    plt.xlabel("Time [s]")
    #plt.ylabel("Phase fraction [-]")
    plt.ylabel("Temperature [degC]")
    plt.ylabel("Plastic strain circumference [-]")
    plt.rcParams.update({'font.size': 30})
    plt.xlim([0,60])
    plt.show()


def DatastreamPlotting(dataname):
    import matplotlib.pyplot as plt

    dataname = "Composition/C"
    filename = "Datastream_Cache.xdmf"
    y = read_results_axis(filename, dataname)
    x = read_results_axis(filename, "nodes")
    plt.plot(x[:, 0], y)
    plt.show()

    names = getnames_results(filename)
    print(names)


    return


if __name__ == "__main__":
    # testing()
    # ResultfileTest()
    # modelling()
    # DatastreamPlotting()
    looping()
    # GUI()
    # DockerTest()
    # ResultPlotting(["Resultfiles/October2024_Ref.xdmf", "Resultfiles/October2024_LPC.xdmf"], "Composition/C")
    # DatastreamPlotting("Composition/C")

    # Result_GUI_show("Resultfiles/October2024.xdmf")
    # Result_GUI_show("Resultfiles/October2024_LPC.xdmf")


    #data = read_input()
    #savedatastream(data["Datastream"]["Savedirect"])