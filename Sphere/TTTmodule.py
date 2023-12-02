import numpy as np

from Sphere.Solvers.Thermocalc import *
from Sphere.Solvers.TTTmodelfit import *
from HelpFile import *

def getTTTcompositions():
    roundingTTT = 2
    data = read_input()
    TTTcompositions = list()
    fullcomposition = dict()
    surfcomp = dict()
    halfcomp = dict()
    for element in data['Material']['Composition'].keys():
        fullcomposition[element] = getaxisvalues("Composition/" + element)
    composition = data['Material']['Composition']
    for element in data['Material']['Composition'].keys():
        surfcomp[element] = round(fullcomposition[element][-1], roundingTTT)
        halfcomp[element] = round((fullcomposition[element][-1] + composition[element]) / 2, roundingTTT)
        tmpcomp = composition.copy()
        tmpcomp[element] = round(fullcomposition[element][-1], roundingTTT)
        TTTcompositions.append(tmpcomp)
    for element in data['Material']['Composition'].keys():
        tmpcomp = surfcomp.copy()
        tmpcomp[element] = round(composition[element], roundingTTT)
        TTTcompositions.append(tmpcomp)
    TTTcompositions.append(composition)
    TTTcompositions.append(halfcomp)
    TTTcompositions.append(surfcomp)
    return TTTcompositions

def runTTTmodule():
    from HelpFile import read_input, checkinput
    if checkinput('TTT'):
        print('Using precalculated TTT simulation')
        return
    print('TTT module')
    TTTcompositions = getTTTcompositions()
    for tmpcomp in TTTcompositions:
        runTTTcalc(tmpcomp)

    return

def runTTTcalc(composition):
    if bool(getTTTdata(composition, "TTTdata")):
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
    if checkinput('ThermoFit'):
        print('Using precalculated phase transformation models')
        return
    print('TTT models module')
    TTTcompositions = getTTTcompositions()
    for tmpcomp in TTTcompositions:
        TTTfit(tmpcomp)
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
    for element in data['Material']['Composition'].keys():
        fullcomposition[element] = readdatastream("Composition/" + element)
    compositions = getTTTcompositions()
    for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
        Z1 = list()
        Z2 =list()
        X = []
        for comp in compositions:
            x = list()
            print(comp)
            TTTdata = getTTTdata(comp, "Modeldata")
            # Get datapoints that aren't Nan
            indx = ~np.isnan(TTTdata[phase][1])
            if not indx.any():
                continue
            if phase == ["Ferrite", "Bainite", "Perlite"]:
                T = TTTdata[phase][0][indx]
                modelpar1 = TTTdata[phase][1][indx]
                modelpar2 = TTTdata[phase][2][indx]
                x.append(list(T))
                m1func = interpolate.splrep(T, modelpar1, s=0)
                m2func = interpolate.splrep(T, modelpar2, s=0)
            else:
                modelpar1 = TTTdata[phase][0][indx] # Ms
                modelpar2 = TTTdata[phase][1][indx] # beta
            # print(T)
            # print(modelpar1)


            # Creating grid for interpolation

            #for element in ["C"]:
            for element in data["Material"]["Composition"].keys():
                tmpx = [comp[element]]*len(T)
                print(tmpx)
                print(x)
                x.append(tmpx)
            x = np.transpose(x)
            print(x)
            input("t")
            if len(X)!=0:
                X = np.vstack([X,x])
            else:
                X = x.copy()
            # TTT data points
            z1 = modelpar1 #interpolate.splev(T, m1func)
            z2 = modelpar2 #interpolate.splev(T, m2func)
            print(z1)
            Z1 = Z1 + list(z1)
            Z2 = Z2 + list(z2)
        interp1 = interpolate.LinearNDInterpolator(X, Z1)
        interp2 = interpolate.LinearNDInterpolator(X, Z2)
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
            saveresult("Modeldata", phase + "/JMAK/T", zdata[0])
            saveresult("Modeldata", phase + "/JMAK/tau",zdata[1])
            saveresult("Modeldata", phase + "/JMAK/n", zdata[2])
        else:
            saveresult("Modeldata", phase + "/T", zdata[0])
            saveresult("Modeldata", phase + "/Ms", zdata[1])
            saveresult("Modeldata", phase + "/beta", zdata[2])
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