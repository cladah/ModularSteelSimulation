from .ModuleStructure_file import CalcModule
from Sphere.HelpFile import saveresult, checkruncondition, getTTTdata, addTTTdata, read_input
from Sphere.Datastream_file import readdatastream, adjustdatastream
from .Solvers.ThermocalcSolver import getTTTcompositions
from .Solvers.TTTmodelfit import JMAKfit, KMfit
import numpy as np


class Transformationmodelmodule(CalcModule):
    def __init__(self):
        super().__init__("Transformationmodels")

    def run(self):
        if not self.runcondition:
            print("Using precalculated " + str(self.module) + " simulation")
            print("Phase transformation model module done\n")
            return

        if self.program == "Python":
            runTTTmodelmodule(self)
            print("Phase transformation model module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

def runTTTmodelmodule(parent):
    print('\nTTT models module')
    if not checkruncondition('ThermoFit'):
        print('Using precalculated phase transformation models')
        return
    TTTcompositions = getTTTcompositions()
    compnr = len(TTTcompositions)
    i = 1
    for tmpcomp in TTTcompositions:
        TTTfit(tmpcomp)
        parent.updateprogress(i / compnr)
        i = i + 1

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

def TTTinterpolatetonodes():
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

    for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
    #for phase in ["Martensite"]:
        Z1 = list()
        Z2 =list()
        X = []
        print("Number of compositions used for interpolations are " + str(len(compositions)))
        for comp in compositions:
            x = list()
            TTTdata = getTTTdata(comp, "Modeldata")

            if phase in ["Ferrite", "Bainite", "Perlite"]:
                gridlen = len(data['Material']['Composition']) + 1 # Adding temp to gridlength
                T = TTTdata[phase][0]
                z1 = TTTdata[phase][1]  # tau
                z2 = TTTdata[phase][2]  # n

                # Checking the values for tau and n if all nan n = 1, tau = -1E12
                if np.isnan(z2).all():
                    z2 = np.nan_to_num(z2, nan=1)
                z2 = np.nan_to_num(z2, nan=float(np.nanmean(z2)))
                z1 = np.nan_to_num(z1, nan=-1E12)
                x.append(list(T))

                # Creating new grid with points for all temperature points
                for element in data["Material"]["Composition"].keys():
                    tmpx = [comp[element]] * len(T)
                    x.append(tmpx)
            else:  # Martensite
                gridlen = len(data['Material']['Composition'])
                z1 = TTTdata[phase][0]  # Ms
                z2 = TTTdata[phase][1]  # beta
                z1 = np.nan_to_num(z1, nan=0.)
                z2 = np.nan_to_num(z2, nan=0.01)
                for element in data["Material"]["Composition"].keys():
                    tmpx = [comp[element]]
                    x.append(tmpx)
            x = np.transpose(x)

            if len(X) != 0:
                X = np.vstack([X, x])
            else:
                X = x.copy()
            # TTT data points
            if phase in ["Ferrite", "Bainite", "Perlite"]:
                Z1 = Z1 + list(z1)
                Z2 = Z2 + list(z2)
            else:
                Z1.append(z1)
                Z2.append(z2)
        # X is all points on the [T, comp1, comp2, comp3, ...] list
        # Z is all points n, tau, Ms, beta values

        # Duplicate check needed if number of elements are reduced
        if phase in ["Ferrite", "Bainite", "Perlite"]:
            dupindx = list()
            for i in range(len(X)):
                for tmpX in X[i+1:]:
                    if (X[i] == tmpX).all():
                        dupindx.append(i)
            X = np.delete(X, dupindx, 0)
            #
            Z1 = np.delete(Z1, dupindx, 0)
            Z2 = np.delete(Z2, dupindx, 0)
            # print(str(len(dupindx)) + " nr of duplicate grid points")
            print(str(len(X)-len(dupindx)) + " nr of grid points used")
        else:
            pass
        tmpX = list(set(np.array(X)[:, 0]))
        tmpX.sort()
        newX = [tmpX]

        # Check len on X grid as T is added for FPB
        for elementnr in range(gridlen-1):
            tmpX = list(set(np.array(X)[:, elementnr+1]))
            tmpX.sort()
            newX.append(tmpX)

        newZ1 = np.array(Z1).reshape(np.shape(np.meshgrid(*newX,indexing='ij')[0]))
        newZ2 = np.array(Z2).reshape(np.shape(np.meshgrid(*newX, indexing='ij')[0]))



        Tgrid = np.linspace(260, 1000, 38)
        grid = list()
        for element in data["Material"]["Composition"].keys():
            if len(grid) == 0:
                grid = [[i] for i in np.array(fullcomposition[element])]
            else:
                tmpgrid = np.array(fullcomposition[element])
                grid = [grid[i] + [tmpgrid[i]] for i in range(len(tmpgrid))]
        z1 = []
        z2 = []

        print("Interpolating modeldata to gridpoints for " + phase)
        #print(np.shape(grid))
        for point in grid:
            if phase in ["Ferrite", "Bainite", "Perlite"]:
                points = [[t] + point for t in Tgrid]
            else:
                points = [point]
            tmpz1 = list()
            tmpz2 = list()
            for p in points:
                indx = points.index(p)
                p = [round(p[i], 1) for i in range(len(p))]
                points[indx] = p
            tmpz1 = tmpz1 + list(interpolate.interpn(newX, newZ1, points))
            tmpz2 = tmpz2 + list(interpolate.interpn(newX, newZ2, points))
            #for p in points:
            #    p = [round(p[i], 1) for i in range(len(p))]
            #    tmpz1 = tmpz1 + list(interpolate.interpn(newX, newZ1, p))
            #    tmpz2 = tmpz2 + list(interpolate.interpn(newX, newZ2, p))
            z1.append(tmpz1)
            z2.append(tmpz2)

        if phase in ["Ferrite", "Bainite", "Perlite"]:
            #z1 = np.nan_to_num(z1,nan=-1E12)
            #z2 = np.nan_to_num(z2, nan=(np.nanmax(z2)+np.nanmin(z2))/2)
            saveresult("Modeldata", phase + "/JMAK/T", Tgrid)
            saveresult("Modeldata", phase + "/JMAK/tau", np.asarray(z1))
            saveresult("Modeldata", phase + "/JMAK/n", np.asarray(z2))
        else:
            saveresult("Modeldata", phase + "/KM/Ms", np.asarray(z1))
            saveresult("Modeldata", phase + "/KM/beta", np.asarray(z2))
        #print(Tgrid)
        #print(np.asarray(z1))

        #input("")

        print(phase + " added to node data")
    print("Modeldata interpolated to nodes")