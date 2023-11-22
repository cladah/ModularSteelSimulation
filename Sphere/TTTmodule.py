from Sphere.Solvers.Thermocalc import *
from Sphere.Solvers.TTTmodelfit import *
from HelpFile import *
from Postprocessing.dataextraction import getaxisvalues


def runTTTmodule():

    from Solvers.Thermocalc import calculateCCT
    import h5py
    import matplotlib.pyplot as plt
    import numpy as np
    from HelpFile import read_input, checkinput
    if checkinput('TTT'):
        print('Using precalculated TTT simulation')
        return
    print('TTT module')
    data = read_input()
    with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
        x = np.array(f.get("CNcurves/Position"))
        a = dict()
        for element in data['Material']['Composition'].keys():
            #plt.plot(x, 100 * np.array(f.get("CNcurves/"+element)))
            a[element] = 100 * np.array(f.get("CNcurves/Elements/" + element)) # Composition curves at all points along x

    composition = data['Material']['Composition']

    runTTTcalc("TTT_center.hdf5", composition)
    composition['C'] = a['C'][-1]
    runTTTcalc("TTT_surface.hdf5", composition)
    return

def runTTTcalc(filename,composition):
    for element in composition:
        saveresult(filename, "Composition/"+element,composition[element])
    Tsteps = np.linspace(270, 1000, 74)
    start, half, finish = calculateMartensite(composition)
    saveresult(filename, "Martensite/start", start)
    saveresult(filename, "Martensite/half", half)
    saveresult(filename, "Martensite/finish", finish)

    start, half, finish = calculatePerlite(Tsteps, composition)
    saveresult(filename, "Perlite/Tsteps", Tsteps)
    saveresult(filename, "Perlite/start", start)
    saveresult(filename, "Perlite/half", half)
    saveresult(filename, "Perlite/finish", finish)

    start, half, finish = calculateBainite(Tsteps, composition)
    saveresult(filename, "Bainite/Tsteps", Tsteps)
    saveresult(filename, "Bainite/start", start)
    saveresult(filename, "Bainite/half", half)
    saveresult(filename, "Bainite/finish", finish)
def runTTTfitmodule():
    if checkinput('ThermoFit'):
        print('Using precalculated phase transformation models')
        return
    print('TTT module')
    files = ["TTT_center.hdf5", "TTT_surface.hdf5"]
    for file in files:
        TTTfit(file)
    TTTinterpolatetonodes(files)
def TTTfit(filename):
    Tlist, n, tau = JMAKfit("Perlite",filename)
    saveresult(filename, "Perlite/JMAK/T", Tlist)
    saveresult(filename, "Perlite/JMAK/n", n)
    saveresult(filename, "Perlite/JMAK/tau", tau)

    b_Tlist, b_n, b_tau = JMAKfit("Bainite",filename)
    saveresult(filename,"Bainite/JMAK/T",b_Tlist)
    saveresult(filename, "Bainite/JMAK/n",b_n)
    saveresult(filename, "Bainite/JMAK/tau",b_tau)

    Ms, beta = KMfit("Martensite",filename)
    saveresult(filename, "Martensite/KM/Ms", Ms)
    saveresult(filename, "Martensite/KM/beta", beta)

def fitspline():
    pass
    #with h5py.File("Resultfiles/TTT.hdf5", "r+") as f:
    #    f.create_dataset('Bainite/JMAK/Tau', data=b_Tlist)
    #   f.create_dataset('Bainite/JMAK/Tau', data=b_n)
    #    f.create_dataset('Bainite/JMAK/Tau', data=b_tau)
def TTTinterpolatetonodes(files):
    from scipy import interpolate
    data = read_input()
    composition = dict()
    for element in data["Material"]["Composition"].keys():
        composition[element] = getaxisvalues(element, 0)
    phases = ["Bainite"]

    for phase in phases:
        phasepm = list()
        for filename in files:
            pm = getJMAK(phase, filename)
            phasepm.append(pm)

    # Testing data
    print(phasepm)
    testT = np.linspace(270, 1000, 74)
    print(interpolate.splev(testT, phasepm[0][1], der=0))
def getJMAK(phase, filename):
    from scipy import interpolate
    T = readresultfile(filename, phase + "/JMAK/T")
    tau = readresultfile(filename, phase + "/JMAK/tau")
    n = readresultfile(filename, phase + "/JMAK/n")
    taufunc = interpolate.splrep(T, tau, s=0)
    nfunc = interpolate.splrep(T, n, s=0)
    return taufunc, nfunc

def getKM():
    pass