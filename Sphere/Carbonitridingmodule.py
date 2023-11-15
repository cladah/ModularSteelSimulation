import numpy as np

from Sphere.Solvers.Thermocalc import *
from HelpFile import *
import h5py
from numpy import interp

def runcarbonitridingmodule():
    if checkinput('Carbonitriding'):
        print('Using precalculated carbnitriding simulation')

        xyz = readdatastream('nodes')
        r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)

        with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
            CN_xyz = np.array(f.get("CNcurves/Position"))
            elements = f.get("CNcurves/Elements").keys()
            for element in elements:
                c_element = np.array(f.get("CNcurves/Elements/" + element))
                elementvalues = interp(r, CN_xyz, c_element)
                adjustdatastream(element, elementvalues, "nodes")
        return
    print('Carbonitriding module')
    data = read_input()
    if data["Programs"]["CNDiffusion"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        CN = TCcarbonitriding(activityenv)
    else:
        raise KeyError(str(data["Programs"]["CNDiffusion"]) + ' not implemented for carbonitriding')

    with h5py.File("Resultfiles/Carbonitriding.hdf5", "w") as f:
        for element in CN[1].keys():
            f.create_dataset("CNcurves/Elements/"+element, data=np.array(CN[1][element]))
        f.create_dataset("CNcurves/Position", data=np.array(CN[0]))



