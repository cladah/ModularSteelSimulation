import numpy as np

from Sphere.Solvers.Thermocalc import *
from HelpFile import *
import h5py
from numpy import interp

class diffusionmodule():
    def __init__(self, run=True):
        data = read_input()
        self.program = data["Programs"]["Carbonitriding"]
        self.run = run
        if checkinput("Carbonitriding"):
            self.run = True
    def reset(self):
        self.run = True
    def run(self):
        print("Carbonitriding module")
        data = read_input()
        if self.run == 0:
            print('Using precalculated carbnitriding simulation')
            for element in data["Material"]["Composition"].keys():
                elementvalues = readdatastreamcache("Composition/" + element)
                adjustdatastream("Composition/" + element, elementvalues, "nodes")
            return
        if self.program == "TC":
            print('Running carbon-nitriding module with ThermoCalc')
            activityenv = TCequalibrium("env")
            composition = TCcarbonitriding(activityenv)

        else:
            raise KeyError(str(data["Programs"]["Carbonitriding"]) + ' not implemented for carbonitriding')

        xyz = readdatastream('nodes')
        r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)
        calc_xyz = np.array(composition[0])
        for element in composition[1].keys():
            calc_value = np.array(composition[1][element])
            nodevalues = interp(r, calc_xyz, calc_value) * 100
            adjustdatastream("Composition/" + element, nodevalues, "nodes")
        print("Carbonitriding modelue done")


def runcarbonitridingmodule():
    print("\nCarbonitriding module")

    xyz = readdatastream('nodes')
    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)

    with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
        CN_xyz = np.array(f.get("CNcurves/Position"))
        elements = f.get("CNcurves/Elements").keys()
        for element in elements:
            c_element = np.array(f.get("CNcurves/Elements/" + element))
            elementvalues = interp(r, CN_xyz, c_element) * 100
            adjustdatastream("Composition/" + element, elementvalues, "nodes")


    data = read_input()
    if checkinput('Carbonitriding'):
        print('Using precalculated carbnitriding simulation')
        for element in data["Material"]["Composition"].keys():
            elementvalues = readdatastreamcache("Composition/" + element)
            adjustdatastream("Composition/" + element, elementvalues, "nodes")
        print("Carbonitriding modelue done")
        return
    if data["Programs"]["Carbonitriding"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        composition = TCcarbonitriding(activityenv)

    else:
        raise KeyError(str(data["Programs"]["Carbonitriding"]) + ' not implemented for carbonitriding')

    xyz = readdatastream('nodes')
    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)
    calc_xyz = np.array(composition[0])
    for element in composition[1].keys():
        calc_value = np.array(composition[1][element])
        nodevalues = interp(r, calc_xyz, calc_value) * 100
        adjustdatastream("Composition/" + element, nodevalues, "nodes")
    print("Carbonitriding modelue done")



def old_runcarbonitridingmodule():
    if checkinput('Carbonitriding'):
        print('Using precalculated carbnitriding simulation')

        xyz = readdatastream('nodes')
        r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)

        with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
            CN_xyz = np.array(f.get("CNcurves/Position"))
            elements = f.get("CNcurves/Elements").keys()
            for element in elements:
                c_element = np.array(f.get("CNcurves/Elements/" + element))
                elementvalues = interp(r, CN_xyz, c_element)*100
                adjustdatastream("Composition/"+element, elementvalues, "nodes")
        return
    print('Carbonitriding module')
    data = read_input()
    if data["Programs"]["Carbonitriding"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        composition = TCcarbonitriding(activityenv)

    else:
        raise KeyError(str(data["Programs"]["Carbonitriding"]) + ' not implemented for carbonitriding')

    # Implement C

    with h5py.File("Resultfiles/Carbonitriding.hdf5", "w") as f:
        for element in composition[1].keys():
            f.create_dataset("CNcurves/Elements/"+element, data=np.array(composition[1][element]))
        f.create_dataset("CNcurves/Position", data=np.array(composition[0]))
    xyz = readdatastream('nodes')
    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)
    with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
        CN_xyz = np.array(f.get("CNcurves/Position"))
        elements = f.get("CNcurves/Elements").keys()
        for element in elements:
            c_element = np.array(f.get("CNcurves/Elements/" + element))
            elementvalues = interp(r, CN_xyz, c_element) * 100
            adjustdatastream("Composition/" + element, elementvalues, "nodes")



