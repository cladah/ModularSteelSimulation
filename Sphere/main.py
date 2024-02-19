import time
import threading
import meshio
import vtk
import pyvista
import numpy as np

from CGUImodule import MainApp
from Datastream_file import createdatastream, createdatastreamcache, removedatastreamcache, savedatastream, \
    readdatastream, readdatastreamcache, getnamesdatastream
from HelpFile import read_input, setupSimulation, createinputcache, change_input, reset_input
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule
from Postprocessing.Postprocess_main import plotcompare

from Testfile import testdatastream, createTTTdiagram_loop
def GUI():
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()

def looping():
    setupSimulation()
    differentin = [["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.6,"Mn": 0.5, "Ni": 1.5,"Mo": 0.3,"Si": 0.2}],
                   ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.0, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
                   ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 0.4,"Mn": 0.5, "Ni": 1.5,"Mo": 0.3,"Si": 0.2}],
                   ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.6,"Mn": 0.5, "Ni": 1.0,"Mo": 0.3,"Si": 0.2}],
                   ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.6,"Mn": 0.5, "Ni": 0.5,"Mo": 0.3,"Si": 0.2}]]
    differentin = [
        ["Material", "Composition", {"C": 0.2, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.3, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.4, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Material", "Composition", {"C": 0.5, "N": 0.025, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3, "Si": 0.2}],
        ["Thermo", "CNtemp", 1000 + 273.15],
        ["Thermo", "CNtemp", 800 + 273.15],
        ["Thermo", "CNtime", 36000/2],
        ["Thermo", "CNtime", 36000/4],
        ["Thermo", "quenchtemp", 223.15], # -50 degC
        ["FEM", "heatflux", {"htc": [400.0,400.0,5600.0,1500.0,1500.0],"T": [0.0,300.0,500.0,650.0,1300.0]}]]  # Double htc

    saveloc = ["Ref.xdmf","C03.xdmf","C04.xdmf","C05.xdmf","CNtemp1000.xdmf", "CNtemp800.xdmf", "CNtime05.xdmf",
               "CNtime025.xdmf", "quenchtemp_n50.xdmf", "heatflux2.xdmf"]

    if len(differentin) != len(saveloc):
        raise KeyError("Wrong looping input")

    i = 0
    for a in differentin:
        print("Loop nr" + str(i+1))
        reset_input()
        change_input(*a)
        data = read_input()
        createdatastreamcache(data["Datastream"]["Cachedirect"])

        # resetdatastream()
        createinputcache()

        modules = list()
        modules.append(Meshingmodule())
        modules.append(Carbonitridingmodule())
        modules.append(TTTdiagrammodule())
        modules.append(Transformationmodelmodule())
        modules.append(Quenchingmodule())

        for currentmodule in modules:
            currentmodule.run()

        removedatastreamcache()
        savedatastream(saveloc[i])
        i = i + 1

def modelling():
    setupSimulation()

    data = read_input()

    createdatastreamcache(data["Datastream"]["Cachedirect"])
    # resetdatastream()
    createinputcache()

    modules = list()
    modules.append(Meshingmodule())
    modules.append(Carbonitridingmodule())
    modules.append(TTTdiagrammodule())
    modules.append(Transformationmodelmodule())
    # modules.append(Quenchingmodule())

    for currentmodule in modules:
        # if currentmodule.modulename() != "Meshing":
        #     tid = threading.Thread(target=run_single_module, args=(currentmodule,))
        #     tid.start()
        #     progressmonitor(tid, currentmodule)  # Making sure thread is done
        # else:
        run_single_module(currentmodule)


    removedatastreamcache()
    savedatastream(data["Datastream"]["Savedirect"])

def run_single_module(module):
    module.run()

def progressmonitor(tid, module):
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)

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

if __name__ == "__main__":
    # modelling()
    # looping()
    GUI()
    # xmdftesting()
    # testdatastream()
    # createTTTdiagram_loop()
    # plotcompare(["Resultfiles/Cr_16.xdmf", "Resultfiles/Cr_10.xdmf", "Resultfiles/Cr_04.xdmf"], "Composition/C", 0)


    # print(np.max(data))
    # data = read_input()
    # createdatastreamcache(data["Datastream"]["Cachedirect"])
    # Meshingmodule().run()
    # Carbonitridingmodule().run()
    # TTTdiagrammodule().run()
    # removedatastreamcache()
    # savedatastream(data["Datastream"]["Savedirect"])


    # Bug in meshio TimeSeriesWriter
    # self.h5_filename = self.filename.stem + ".h5"
    # self.h5_filename = self.filename.with_suffix(".h5")