import time

import meshio

from CGUImodule import MainApp
from Datastream_file import createdatastream, createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input, setupSimulation, createinputcache
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


def modelling():
    setupSimulation()

    data = read_input()
    import threading
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

if __name__ == "__main__":
    modelling()
    #Testfile.read_data_from_xdmf("Resultfiles/Test.xdmf", 0)
    #Testfile.add_data_to_xdmf("Resultfiles/Datastream.xdmf", [], 0)

    #modelling()
    # data = read_input()
    # createdatastreamcache(data["Datastream"]["Cachedirect"])
    # Meshingmodule().run()
    # Carbonitridingmodule().run()
    # TTTdiagrammodule().run()
    # removedatastreamcache()
    # savedatastream(data["Datastream"]["Savedirect"])