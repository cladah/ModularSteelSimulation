import time

import numpy as np

from Sphere.Solvers.ThermocalcSolver import *
from HelpFile import *
import h5py
from numpy import interp

class diffusionmodule():
    def __init__(self, run=True):
        data = read_input()
        self.program = data["Programs"]["Carbonitriding"]
        self.run = run
        if checkruncondition("Carbonitriding"):
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
        print("Carbonitriding module done")


def runcarbonitridingmodule(parent):
    print("Carbonitriding module")
    parent.updateprogress(0.1)

    data = read_input()
    if not checkruncondition('Carbonitriding'):
        print('Using precalculated carbnitriding simulation')
        for element in data["Material"]["Composition"].keys():
            elementvalues = readdatastreamcache("Composition/" + element)
            adjustdatastream("Composition/" + element, elementvalues, "nodes")
        print("Carbonitriding module done")
        parent.updateprogress(1.0)
        return

    if data["Programs"]["Carbonitriding"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        print("Avtivity of atmosphere calculated")
        parent.updateprogress(0.2)
        composition = TCcarbonitriding(activityenv)
        parent.updateprogress(0.9)
    else:
        raise KeyError(str(data["Programs"]["Carbonitriding"]) + ' not implemented for carbonitriding')

    xyz = readdatastream('nodes')
    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)
    calc_xyz = np.array(composition[0])
    for element in composition[1].keys():
        calc_value = np.array(composition[1][element])
        nodevalues = interp(r, calc_xyz, calc_value) * 100
        adjustdatastream("Composition/" + element, nodevalues, "nodes")
    parent.updateprogress(1.0)
    print("Carbonitriding module done")