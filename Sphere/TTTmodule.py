import numpy as np

from Sphere.Solvers.Thermocalc import *
from Sphere.Solvers.TTTmodelfit import *
from HelpFile import *

def getTTTcompositions():
    roundingTTT = 1
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
        if tmpcomp not in TTTcompositions:
            TTTcompositions.append(tmpcomp)
    if composition not in TTTcompositions:
        TTTcompositions.append(composition)
    if halfcomp not in TTTcompositions:
        TTTcompositions.append(halfcomp)
    if surfcomp not in TTTcompositions:
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
    # Take away values that are 0 EXCEPT N AND C!
    keys = composition.keys()
    for key in keys:
        if key in ["C", "N"]:
            continue
        if composition[key] == 0.:
            del composition[key]
    #
    if bool(getTTTdata(composition, "TTTdata")):
        print("TTTdata exists in database for " + str(composition))
        return
    print("Running TTT calculation for " + str(composition))
    Tsteps = np.linspace(300, 1000, 71)  # 71
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
            start = [[start, start], [0.1,1E12]]
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
    newcomp = list()
    for comp in compositions:
        if comp not in newcomp:
            newcomp.append(comp)
    compositions = newcomp.copy()

    for phase in ["Ferrite", "Bainite", "Perlite"]:
        Z1 = list()
        Z2 =list()
        X = []
        for comp in compositions:
            x = list()
            TTTdata = getTTTdata(comp, "Modeldata")
            # Get datapoints that aren't Nan
            indx = ~np.isnan(TTTdata[phase][1])
            if not indx.any():
                continue

            if phase in ["Ferrite", "Bainite", "Perlite"]:
                T = TTTdata[phase][0][indx]
                z1 = TTTdata[phase][1][indx]
                z2 = TTTdata[phase][2][indx]
                x.append(list(T))
                for element in data["Material"]["Composition"].keys():
                    tmpx = [comp[element]] * len(T)
                    x.append(tmpx)
            else:
                z1 = TTTdata[phase][0][indx] # Ms
                z2 = TTTdata[phase][1][indx] # beta
                for element in data["Material"]["Composition"].keys():
                    tmpx = [comp[element]]
                    x.append(tmpx)
            x = np.transpose(x)

            if len(X)!=0:
                X = np.vstack([X, x])
            else:
                X = x.copy()
            # TTT data points
            Z1 = Z1 + list(z1)
            Z2 = Z2 + list(z2)
        # X is all points on the [T, comp1, comp2, comp3, ...] list
        # Z is all points n, tau, Ms, beta values

        # Reducing the amount of elements
        X = X[:,0:3]
        # Duplicate check needed if number of elements are reduced
        if phase in ["Ferrite", "Bainite", "Perlite"]:
            dupindx = list()
            for i in range(len(X)):
                for tmpX in X[i+1:]:
                    if (X[i]==tmpX).all():
                        dupindx.append(i)
            X = np.delete(X, dupindx, 0)
            #
            Z1 = np.delete(Z1, dupindx, 0)
            Z2 = np.delete(Z2, dupindx, 0)
            #print(str(len(dupindx)) + " nr of duplicate compositions")

        interp1 = interpolate.LinearNDInterpolator(X, Z1)
        interp2 = interpolate.LinearNDInterpolator(X, Z2)

        Tgrid = np.linspace(300, 1000, 71)
        #print(Tgrid)
        grid = list()
        for element in data["Material"]["Composition"].keys():
            if len(grid) == 0:
                grid = [[i] for i in np.array(fullcomposition[element])]
            else:
                tmpgrid = np.array(fullcomposition[element])
                grid = [grid[i] + [tmpgrid[i]] for i in range(len(tmpgrid))]
        z1 = list()
        z2 = list()


        #print(grid[:, 0])
        #print([grid[i][0:2] for i in range(len(grid))])
        grid = [grid[i][0:2] for i in range(len(grid))]
        for point in grid:
            point = [[t] + point for t in Tgrid]

            z1.append(interp1(point))
            z2.append(interp2(point))
        print("Results has the shape")
        print(np.shape(Tgrid))
        print(np.shape(z1))
        print(np.shape(z2))
        if phase in ["Ferrite", "Bainite", "Perlite"]:
            z1 = np.nan_to_num(z1,nan=-1E12)
            z2 = np.nan_to_num(z2, nan=(np.nanmax(z2)+np.nanmin(z2))/2)
            saveresult("Modeldata", phase + "/JMAK/T", Tgrid)
            saveresult("Modeldata", phase + "/JMAK/tau", np.asarray(z1))
            saveresult("Modeldata", phase + "/JMAK/n", np.asarray(z2))
        else:
            saveresult("Modeldata", phase + "/Ms", np.asarray(z1))
            saveresult("Modeldata", phase + "/beta", np.asarray(z2))
        print(phase + " added to node data")
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