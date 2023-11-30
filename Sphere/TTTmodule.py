import numpy as np

from Sphere.Solvers.Thermocalc import *
from Sphere.Solvers.TTTmodelfit import *
from HelpFile import *


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
        fullcomposition[element] = getaxisvalues("Composition/" + element)
    composition = data['Material']['Composition']
    surfcomp = dict()
    halfcomp = dict()
    TTTcompositions = list()
    for element in data['Material']['Composition'].keys():
        surfcomp[element] = round(fullcomposition[element][-1], 2)
        halfcomp[element] = round((fullcomposition[element][-1]+composition[element])/2, 2)
        tmpcomp = composition.copy()
        tmpcomp[element] = round(fullcomposition[element][-1], 2)
        TTTcompositions.append(tmpcomp)

    runTTTcalc(composition)
    runTTTcalc(halfcomp)
    runTTTcalc(surfcomp)
    for tmpcomp in TTTcompositions:
        runTTTcalc(tmpcomp)

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
def runTTTmodelmodule():
    import h5py
    if checkinput('ThermoFit'):
        print('Using precalculated phase transformation models')
        return
    print('TTT models module')
    data = read_input()
    composition = data['Material']['Composition']
    surfcomp = dict()
    for element in composition.keys():
        surfcomp[element] = round(getaxisvalues("Composition/" + element)[-1], 2)
    TTTfit(composition)
    TTTfit(surfcomp)
    print("Models fitted to data")
    TTTinterpolatetonodes()

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
    import matplotlib.pyplot as plt
    from scipy import interpolate
    data = read_input()
    fullcomposition = dict()
    fullcomposition = dict()
    core = dict()
    surface = dict()
    complist = list()
    phase = "Bainite"
    for element in data["Material"]["Composition"].keys():
        compaxis = getaxisvalues("Composition/"+element)
        core[element] = compaxis[0]
        surface[element] = compaxis[-1]
        fullcomposition[element] = readdatastream("Composition/"+element)
    coremodel = getTTTdata(core, "Modeldata")
    surfacemodel = getTTTdata(surface, "Modeldata")
    for phase in ["Bainite", "Perlite", "Martensite"]:
        if phase in ["Ferrite","Bainite","Perlite"]:
            indx = ~np.isnan(coremodel[phase][1])
            T = coremodel[phase][0][indx]
            tau = coremodel[phase][1][indx]
            n = coremodel[phase][2][indx]
            taufunc = interpolate.splrep(T, tau, s=0)
            nfunc = interpolate.splrep(T, n, s=0)
            indx = ~np.isnan(surfacemodel[phase][1])
            T2 = surfacemodel[phase][0][indx]
            tau2 = surfacemodel[phase][1][indx]
            n2 = surfacemodel[phase][2][indx]
            taufunc2 = interpolate.splrep(T2, tau2, s=0)
            nfunc2 = interpolate.splrep(T2, n2, s=0)
            complist =[list(core.values()),list(surface.values())]

            # Creating grid for JMAK interpolation
            x = list()
            grid = list()
            x.append(np.append(T,T2))
            for element in ["C"]:
            #for element in data["Material"]["Composition"].keys():
                tmpx = [core[element]]*len(T) + [surface[element]]*len(T2)
                x.append(tmpx)
                grid.append((min(tmpx) + max(tmpx))/2)
            x = np.transpose(x)

            # TTT data points
            z = np.append(interpolate.splev(T, taufunc), interpolate.splev(T2, taufunc2))
            interp1 = interpolate.LinearNDInterpolator(x, z)
            z = np.append(interpolate.splev(T, nfunc), interpolate.splev(T2, nfunc2))
            interp2 = interpolate.LinearNDInterpolator(x, z)
            intermodels = [interp1, interp2]

        Tgrid = T
        grid = np.array(fullcomposition["C"])
        zdata = list()
        for point in grid:
            point = [[t,point] for t in Tgrid]
            for inter in intermodels:
                Z = inter(point)
            zdata.append(Z)
        zdata = np.asarray(zdata)
        if phase in ["Ferrite", "Bainite", "Perlite"]:
            saveresult("Modeldata", phase + "/JMAK/tau",zdata)
            saveresult("Modeldata", phase + "/JMAK/n", zdata)
        else:
            saveresult("Modeldata", phase + "/Ms", zdata)
            saveresult("Modeldata", phase + "/beta", zdata)
    print("Modeldata interpolated to nodes")
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