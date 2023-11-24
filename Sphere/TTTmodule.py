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
    fullcomposition = dict()
    for element in data['Material']['Composition'].keys():
        fullcomposition[element] = getaxisvalues("N")[1]
    composition = data['Material']['Composition']
    surfcomp = dict()
    for element in data['Material']['Composition'].keys():
        surfcomp[element] = round(fullcomposition[element][-1], 2)
    runTTTcalc(composition)
    runTTTcalc(surfcomp)
    return

def runTTTcalc(composition):
    if bool(getTTTdata(composition,"TTTdata")):
        print("TTTdata exists in database for " + str(composition))
        return
    print("Running TTT calculation for " + str(composition))
    Tsteps = np.linspace(300, 1000, 21)  # 74
    phases = ["Ferrite", "Bainite", "Perlite","Martensite"]
    TTTdata = dict()
    for ph in phases:
        if ph == "Ferrite":
            start, half, finish = calculateFerrite(Tsteps, composition)
            start = [Tsteps, start]
            half = [Tsteps, half]
            finish = [Tsteps, finish]
        elif ph == "Perlite":
            start, half, finish = calculatePerlite(Tsteps, composition)
            start = [Tsteps, start]
            half = [Tsteps, half]
            finish = [Tsteps, finish]
        elif ph == "Bainite":
            start, half, finish = calculateBainite(Tsteps, composition)
            start = [Tsteps, start]
            half = [Tsteps, half]
            finish = [Tsteps, finish]
        elif ph == "Martensite":
            start, half, finish = calculateMartensite(composition)
            start = [[start,start],[0.1,1E12]]
            half = [[half, half], [0.1, 1E12]]
            finish = [[finish, finish], [0.1, 1E12]]
        else:
            raise KeyError("Phase error in TTT module, check input phases")
        phase = dict()
        phase["start"] = start
        phase["half"] = half
        phase["finish"] = finish
        TTTdata[ph] = phase
    addTTTdata(composition, TTTdata, "TTTdata")
def runTTTfitmodule():
    import h5py
    if checkinput('ThermoFit'):
        print('Using precalculated phase transformation models')
        return
    print('TTT fitting module')
    data = read_input()
    composition = data['Material']['Composition']
    surfcomp = dict()
    for element in composition.keys():
        surfcomp[element] = round(getaxisvalues(element, 0)[-1],2)
    TTTfit(composition)
    TTTfit(surfcomp)
    print("models fitted to data")
    #TTTinterpolatetonodes()
def TTTfit(composition):
    phases = ["Ferrite","Perlite","Bainite","Martensite"]
    modeldata = dict()
    for phase in phases:
        if phase in ["Ferrite", "Perlite", "Bainite"]:
            T, tau, n = JMAKfit(composition, phase)
            modeldata[phase] = [T, tau, n]
        elif phase == "Martensite":
            Ms, beta = KMfit(composition, phase)
            modeldata[phase] = [Ms, beta]
    addTTTdata(composition, modeldata, "Modeldata")
def fitspline():
    pass
    #with h5py.File("Resultfiles/TTT.hdf5", "r+") as f:
    #    f.create_dataset('Bainite/JMAK/Tau', data=b_Tlist)
    #   f.create_dataset('Bainite/JMAK/Tau', data=b_n)
    #    f.create_dataset('Bainite/JMAK/Tau', data=b_tau)
def TTTinterpolatetonodes():
    data = read_input()
    for element in data["Material"]["Composition"].keys():
        readdatastream("Composition/"+element)


    #getTTTdata(composition)
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
    #print(phasepm)
    testT = np.linspace(270, 1000, 74)
    print(interpolate.splev(testT, phasepm[0][1], der=0))
def getJMAK(composition, phase):
    from scipy import interpolate
    T = readresultfile(filename, phase + "/JMAK/T")

    tau = readresultfile(filename, phase + "/JMAK/tau")

    n = readresultfile(filename, phase + "/JMAK/n")
    # Taking out Nan values
    xpoints = -tau * (-np.log(0.98)) ** (1 / n)
    indx = xpoints < 1E12 * np.logical_not(np.isnan(tau))
    T = T[indx]
    tau = tau[indx]
    n = n[indx]
    taufunc = interpolate.splrep(T, tau, s=0)
    nfunc = interpolate.splrep(T, n, s=0)
    return taufunc, nfunc

def getKM():
    pass