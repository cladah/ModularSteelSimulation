from Sphere.Solvers.Thermocalc import *
from HelpFile import *
import h5py
from numpy import interp

def runcarbonitriding():
    if checkinput('Carbonitriding'):
        print('Using precalculated carbnitriding simulation')
        return
    print('Carbonitriding module')
    data = read_input()
    if data["Programs"]["CNDiffusion"] == "TC":
        print('Running carbon-nitriding module with ThermoCalc')
        activityenv = TCequalibrium("env")
        CN = TCcarbonitriding(activityenv)
    else:
        raise KeyError(str(data["Programs"]["CNDiffusion"]) + ' not implemented for carbonitriding')

    xyz = readdatastream('nodes')
    r = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2 + xyz[:, 2] ** 2)

    with h5py.File("Resultfiles/Carbonitriding.hdf5", "w") as f:
        for element in CN[1].keys():
            f.create_dataset("CNcurves/"+element, data=np.array(CN[1][element]))
        f.create_dataset("CNcurves/Position", data=np.array(CN[0]))

    for element in CN[1].keys():
        elementvalues = interp(r,np.array(CN[0]),np.array(CN[1][element]))
        adjustdatastream(element, elementvalues, "nodes")

