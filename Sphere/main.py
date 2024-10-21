import time
import meshio
import numpy as np
from matplotlib import pyplot as plt

from CGUImodule import MainApp
from Result_GUI import Result_MainApp, Compare_MainApp
from Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input, setupSimulation, createinputcache, change_input, reset_input, analyseTTTdatabase
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
        savedatastream(data["Datastream"]["Savedirect"])

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


def ResultPlotting(filenames, dataname, point=[0., 0.], tid=-1):
    import matplotlib.pyplot as plt
    for filename in filenames:
        names = getnames_results(filename)
        print("Results in " + filename)
        print(names)
        tmpdata = read_results_axis(filename, dataname, tid)
        xyz = read_results_axis(filename, "nodes")
        plt.plot(xyz[:, 0], tmpdata)
    plt.legend(filenames)
    plt.title(dataname)
    #plt.xlabel("Time [s]")
    #plt.ylabel("Temperature [degC]")
    plt.rcParams.update({'font.size': 30})
    plt.show()
    print("Done")


def checkDB():
    analyseTTTdatabase()


def DatastreamPlotting(dataname):
    import matplotlib.pyplot as plt
    filename = "Datastream.xdmf"
    print(getnames_results("Datastream.xdmf"))
    y = read_results_axis(filename, dataname)
    x = read_results_axis(filename, "nodes")
    plt.plot(x[:, 0], y)
    plt.show()

def TCtest():
    version = '2022b'

    print('Testing TC-Python version: ' + version)
    print(
        'Please make sure that the variable "version" above, matches the release that you want to test, if not change it and re-run this script.')

    # below this line, nothing needs to be manually updated.

    import sys
    print('')
    print('Python version: (should be at least 3.5 and can NOT be older than 3.0)')
    print(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))
    if sys.version_info[0] < 3 or sys.version_info[1] < 5:
        print('Wrong version of Python !!!!!')

    print('')
    print(
        'Python executable path: (gives a hint about the used virtual / conda environment, in case of Anaconda the corresponding \n'
        'environment name can be found by running `conda env list` on the Anaconda command prompt, '
        'TC-Python must be installed into \nEACH separate environment used!)')
    print(sys.executable)

    import os
    print('')
    print(
        'Thermo-Calc ' + version + ' installation directory: (must be a valid path to a complete installation of ' + version + ')')
    tc_env_variable = 'TC' + version[2:].upper() + '_HOME'
    try:
        print(os.environ[tc_env_variable])
    except:
        print('No Thermo-calc environment variable for ' + version + ' was found. (' + tc_env_variable + ')')

    print('')
    print('Url of license server: (if license server is NO-NET, you need a local license file)')
    try:
        print(os.environ['LSHOST'])
    except:
        print('No Thermo-calc license server url was found. (LSHOST)')

    print('')
    print('Path to local license file: (only necessary if not using license server)')
    try:
        print(os.environ['LSERVRC'])
    except:
        print('No path to local license file was found. (LSERVRC)')

    import tc_python
    numerical_version = version[:-1]
    if version[-1] == 'a':
        numerical_version += '.1.*'
    elif version[-1] == 'b':
        numerical_version += '.2.*'
    print('')
    print('TC-Python version: (needs to be ' + numerical_version + ')')
    print(tc_python.__version__)

    with tc_python.TCPython() as session:
        print('')
        print(
            'Lists the databases: (should be a complete list of the installed databases that you have license for or do not require license)')
        print(session.get_databases())

if __name__ == "__main__":
    # testing()
    # ResultfileTest()
    # TCtest()
    # modelling()
    # DatastreamPlotting("Composition/C")
    # print(getnames_results("Resultfiles/October2024_ref.xdmf"))
    checkDB()

    # looping()
    # GUI()
    # DockerTest()

    #dataname = "Composition/C"
    #ResultPlotting(["Resultfiles/October2024_900C.xdmf", "Resultfiles/October2024_Ref.xdmf",
    #                "Datastream.xdmf", "Resultfiles/October2024_LPC.xdmf","Resultfiles/October2024_LPC_4h.xdmf"], dataname)

    # DatastreamPlotting("Composition/C")

    # Result_GUI_show("Resultfiles/October2024_900C.xdmf")
    # Result_GUI_show("Resultfiles/October2024_LPC.xdmf")


    #data = read_input()
    #savedatastream(data["Datastream"]["Savedirect"])