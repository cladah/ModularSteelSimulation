import time
import threading
import meshio
import vtk
import pyvista
import numpy as np

from CGUImodule import MainApp
from Datastream_file import createdatastream, createdatastreamcache, removedatastreamcache, savedatastream, \
    readdatastream, readdatastreamcache, getnamesdatastream
from HelpFile import read_input, setupSimulation, createinputcache, change_input
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule

import Testfile

def GUI():
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()

def looping():
    setupSimulation()
    differentin = [["Material", "Composition", {"C": 0.2,"N": 0.025,"Cr": 1.6,"Mn": 0.5,"Ni": 1.5,"Mo": 0.3,"Si": 0.2}],
                   ["Material", "Composition", {"C": 0.5,"N": 0.025,"Cr": 1.6,"Mn": 0.5,"Ni": 1.5,"Mo": 0.3,"Si": 0.2}],
                   ["Material", "Composition", {"C": 1.0,"N": 0.025,"Cr": 1.6,"Mn": 0.5,"Ni": 1.5,"Mo": 0.3,"Si": 0.2}]]
    saveloc = ["Test1.xmdf", "Test2.xmdf", "Test3.xmdf"]
    i = 0
    for a in differentin:
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
            if currentmodule.modulename() != "Meshing":
                tid = threading.Thread(target=run_single_module, args=(currentmodule,))
                tid.start()
                progressmonitor(tid, currentmodule)  # Making sure thread is done
            else:
                run_single_module(currentmodule)

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
    modules.append(Quenchingmodule())

    for currentmodule in modules:
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=run_single_module, args=(currentmodule,))
            tid.start()
            progressmonitor(tid, currentmodule)  # Making sure thread is done
        else:
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
    data = read_input()
    createdatastreamcache(data["Datastream"]["Cachedirect"])

    with meshio.xdmf.TimeSeriesReader("Datastream_Cache.xdmf") as reader:
        points, cells = reader.read_points_cells()
        pd_list, cd_list, t_list = list(), list(), list()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            t_list.append(t)
            pd_list.append(point_data)
            cd_list.append(cell_data)
        print("Points")
        print(cells[0].data)
        print(cell_data)
        print("Data list")
        print(pd_list[1].keys())
        #print(readdatastream("Austenite", time=-1))

if __name__ == "__main__":
    #Testfile.read_data_from_xdmf("Resultfiles/230124_2.xdmf", 0)
    #Testfile.read_data_from_xdmf("Resultfiles/Test.xdmf", 0)
    #Testfile.add_data_to_xdmf("Resultfiles/Datastream.xdmf", [], 0)

    modelling()
    # GUI()
    # xmdftesting()



    #print(np.max(data))
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