import time

from CGUImodule import MainApp
from Datastream_file import createdatastream, createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule


def GUI():
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()


def modelling():
    data = read_input()
    import threading
    createdatastreamcache(data["Datastream"]["Cachedirect"])

    modules = list()
    modules.append(Meshingmodule())
    modules.append(Carbonitridingmodule())
    modules.append(TTTdiagrammodule())
    #modules.append(CalcModule("TTTmodeling"))
    #modules.append(CalcModule("Quenching"))

    for currentmodule in modules:
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=run_single_module, args=(currentmodule,))
            tid.start()
            progressmonitor(tid, currentmodule)
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
    # data = read_input()
    # createdatastreamcache(data["Datastream"]["Cachedirect"])
    # Meshingmodule().run()
    # Carbonitridingmodule().run()
    # TTTdiagrammodule().run()
    # removedatastreamcache()
    # savedatastream(data["Datastream"]["Savedirect"])