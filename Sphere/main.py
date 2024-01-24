import time

from CGUImodule import MainApp
from HelpFile import read_input, createdatastreamcache, removedatastreamcache, savedatastream
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule
from Sphere.Modulefiles.Tempering_file import Temperingmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
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
    modules.append(Meshingmodule)
    modules.append(Carbonitridingmodule)
    #modules.append(CalcModule("Meshing"))
    #modules.append(CalcModule("Carbonitriding"))
    #modules.append(CalcModule("TTT"))
    #modules.append(CalcModule("TTTmodeling"))
    #modules.append(CalcModule("Quenching"))

    for currentmodule in modules:
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=run_single_module, args=(currentmodule,))
            tid.start()
            progressmonitor(tid, currentmodule)
        else:
            run_single_module(currentmodule)
    raise KeyError
    input("")
    removedatastreamcache()
    savedatastream(data["Datastream"]["Savedirect"])

def run_single_module(module):
    module.run()

def progressmonitor(tid, module):
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)

if __name__ == "__main__":
    data = read_input()
    createdatastreamcache(data["Datastream"]["Cachedirect"])
    Meshingmodule().run()
    Carbonitridingmodule().run()
    TTTdiagrammodule().run()
    removedatastreamcache()
    savedatastream(data["Datastream"]["Savedirect"])