from .ModuleStructure_file_new import CalcModule
from Framework.HelpFile import saveresult, checkruncondition, getTTTdb, addTTTdb, read_input, read_geninput, checkTTTdb
from Framework.Datastream_file import readdatastream, adjustdatastream, readdatastreamcache
from .Solvers.ThermocalcSolver import getTTTcompositions
from .Solvers.TTTmodelfit import JMAKfit, KMfit
import numpy as np


class Transformationmodelmodule(CalcModule):
    def __init__(self, infile, modulenr):
        infile = "Inputs/" + infile + ".json"
        super().__init__("TransformMod", infile, modulenr)

    def run(self):
        models  = {"JMAK":"1-exp(-\u03C4t^n)","KM":"1-exp(-\u03B2(Ms-T))"}
        outstr = ["\n---------------------------------------------------------------------\n"
                  "Transformation model function module: " + self.inputfile + "\n"]
        for ph in self.minput["Phases"]:
            tmpmod = models[self.minput["Model"][ph]]
            outstr.append(f"{ph} modeled with {tmpmod}")
        outstr.append("\n---------------------------------------------------------------------\n\n")


        for line in outstr:
            self.writeoutput(line)
            print(line)

        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")
            precalcinp = ["JMAK_tau_Ferrite", "JMAK_tau_Pearlite", "JMAK_tau_Bainite", "JMAK_n_Ferrite",
                          "JMAK_n_Pearlite",
                          "JMAK_n_Bainite", "KM_Ms_Martensite", "KM_b_Martensite"]
            for pre in precalcinp:
                values = readdatastreamcache(pre)
                adjustdatastream({pre: values}, datapos="nodes")
            print("Phase transformation model module done\n")
            return

        if self.program == "Python":
            runTTTmodelmodule(self)
            print("Phase transformation model module done\n")
        else:
            raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")

def runTTTmodelmodule(parent):
    print('\nTTT models module')
    TTTcompositions = getTTTcompositions()
    compnr = len(TTTcompositions)
    i = 1
    for tmpcomp in TTTcompositions:
        TTTfit(tmpcomp, parent.minput)
        parent.updateprogress(i / compnr)
        i = i + 1

    print("Models fitted to data")
    TTTpolyfit(parent.minput)

    # TTTinterpolatetonodes()

def TTTfit(composition, minput):
    phases = ["Ferrite","Pearlite","Bainite","Martensite"]
    phases = list(minput["Phases"])

    model_input = {k: minput[k] for k in ["GrainSize"]}

    for phase in phases:
        growth_type = minput["GrowthMode"][phase]
        model_type = minput["Model"][phase]
        if not checkTTTdb(composition, model_input, f"{phase}_{growth_type}"):
            next

        raw_data = getTTTdb(composition, model_input, f"{phase}_{growth_type}")

        if phase in ["Ferrite", "Pearlite", "Bainite"]:
            T, tau, n = JMAKfit(raw_data, minput)
            modeldata = [T, tau, n]
        elif phase == "Martensite":
            Ms, beta = KMfit(raw_data, minput)
            modeldata = [Ms, beta]
        else:
            raise KeyError("Phase entered into input not implemented in TTTfit")
        addTTTdb(modeldata, composition, model_input, f"{phase}_{growth_type}_{model_type}")

def TTTinterpolatetonodes():
    from scipy import interpolate
    data = read_geninput()

    fullcomposition = dict()
    for element in data['Material']['Composition'].keys():
        fullcomposition[element] = readdatastream("Composition/" + element)
    compositions = getTTTcompositions()
    newcomp = list()
    for comp in compositions:
        if comp not in newcomp:
            newcomp.append(comp)
    compositions = newcomp.copy()

    for phase in ["Ferrite", "Bainite", "Pearlite", "Martensite"]:
    # for phase in ["Martensite"]:
        Z1 = list()
        Z2 =list()
        X = []
        print("Number of compositions used for interpolations are " + str(len(compositions)))
        for comp in compositions:
            x = list()
            TTTdata = getTTTdb(comp, model_input, f"{phase}_{growth_type}")
            if phase in ["Ferrite", "Bainite", "Pearlite"]:
                gridlen = len(data['Material']['Composition']) + 1 # Adding temp to gridlength
                T = TTTdata[phase][0]
                z1 = TTTdata[phase][1]  # tau
                z2 = TTTdata[phase][2]  # n


                # Checking the values for tau and n if all nan n = 3, tau = 1E12
                if np.isnan(z2).all():
                    z2 = np.nan_to_num(z2, nan=3.)
                #z2 = np.nan_to_num(z2, nan=float(np.nanmean(z2)))
                z2 = np.nan_to_num(z2, nan=3.)

                z2[z2 > 0] = 3.

                z1 = np.nan_to_num(z1, nan=1E12)
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
            if phase in ["Ferrite", "Bainite", "Pearlite"]:
                Z1 = Z1 + list(z1)
                Z2 = Z2 + list(z2)
            else:
                Z1.append(z1)
                Z2.append(z2)
        # X is all points on the [T, comp1, comp2, comp3, ...] list
        # Z is all points n, tau, Ms, beta values

        # Duplicate check needed if number of elements are reduced
        if phase in ["Ferrite", "Bainite", "Pearlite"]:
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

        # print(newX)
        # print(Z1[1])

        newZ1 = np.transpose(np.array(Z1).reshape(np.shape(np.transpose(np.meshgrid(*newX, indexing='ij')[0]))))
        newZ2 = np.transpose(np.array(Z2).reshape(np.shape(np.transpose(np.meshgrid(*newX, indexing='ij')[0]))))
        # print(newZ1[1,0,0,0,0,0,0])
        # input("")

        # import matplotlib.pyplot as plt
        # indx = [i for i, v in enumerate(Z1[0:41]) if v < 1E12]
        # indx = [np.min(indx)-1]+indx + [np.max(indx)+1]
        # x = X[0:41, 0][indx]
        # z1 = np.polyfit(x, np.log(Z1[0:41][indx]), 2)
        # p1 = np.poly1d(z1)
        # z2 = np.polyfit(x, np.log(Z2[0:41][indx]), 2)
        # p2 = np.poly1d(z2)
        # print(Z1[0:41][indx])
        # print(X[0:41, 0][indx])
        # #x = np.linspace(0,1000,40)
        # y = np.array(np.exp(p1(x)))*(-np.log(0.02))**np.array(np.exp(p2(x)))
        # Y = np.array(Z1[0:41])*(-np.log(0.02))**np.array(Z2[0:41])
        #
        # plt.plot(y, x)
        # plt.plot(Y, X[0:41, 0])
        # plt.xscale("log")
        # plt.show()
        # input("STOP")

        #Tgrid = np.linspace(0, 975, 40) + 273.15
        Tgrid = np.linspace(0, 975, 20) + 273.15
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
        newgrid = list()
        for point in grid:
            if phase in ["Ferrite", "Bainite", "Pearlite"]:

                points = [[t] + point for t in Tgrid]
                tmppoints = list()
                for tmppoint in points:
                    tmppoint = [tmppoint[0]] + [round(tmppoint[k], 4) for k in range(1, 3)] + [round(tmppoint[k], 1) for k in
                                                                       range(3, len(tmppoint))]
                    tmppoints.append(tmppoint)
                newgrid.append(tmppoints)
            else: # Martensite
                tmppoint = [round(point[k], 4) for k in range(0, 2)] + [round(point[k], 1) for k in
                                                                       range(2, len(point))]
                newgrid.append([tmppoint])


        # print(np.shape(newgrid))
        # print(newX)
        # print(np.shape(newZ1))
        z1 = list(interpolate.interpn(newX, newZ1, newgrid))
        z2 = list(interpolate.interpn(newX, newZ2, newgrid))
        #print(z2)
        if phase in ["Ferrite", "Bainite", "Pearlite"]:
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

        print(phase + " added to Modeldata file")
    print("Modeldata interpolated to nodes")


import numpy as np
from scipy import interpolate
from pathlib import Path


def TTTpolyfit(minput):
    data = read_geninput()
    model_input = {k: minput[k] for k in ["GrainSize"]}

    comp_elements = list(data['Material']['Composition'].keys())

    # 1. Fetch Composition for all nodes (Vectorized)
    full_composition = {el: readdatastream(f"Composition_{el}") for el in comp_elements}
    # 2. Get unique compositions from database
    raw_compositions = getTTTcompositions()
    # Use a set of tuples to find unique compositions while preserving element maps
    unique_compositions = []
    seen = set()
    for c in raw_compositions:
        c_tuple = tuple(c.items())
        if c_tuple not in seen:
            unique_compositions.append(c)
            seen.add(c_tuple)

    print(f"Number of compositions for interpolation: {len(unique_compositions)}")

    # 3. Process each phase
    phases = list(minput["Phases"])
    for phase in phases:
        growth_type = minput["GrowthMode"][phase]
        model_type = minput["Model"][phase]

        Z1_list, Z2_list, X_list = [], [], []

        for comp in unique_compositions:
            ttt_data = getTTTdb(comp, model_input, f"{phase}_{growth_type}_{model_type}")
            # Extract coordinates for interpolation grid
            X_list.append([comp[el] for el in comp_elements])

            if phase != "Martensite":
                # JMAK Logic (Ferrite, Bainite, Pearlite)
                T, tau, n = ttt_data[0:3]

                # Robust NaN handling
                taumax = 1e8
                n = np.nan_to_num(n, nan=3.0)
                tau = np.nan_to_num(tau, nan=taumax)
                tau = np.clip(tau, None, taumax)

                indices = np.where(tau < taumax)[0]

                if len(indices) == 0:
                    # Default flat profile if no valid TTT data
                    z1 = np.append(np.zeros(5), taumax)
                    z2 = np.append(np.zeros(5), 3.0)
                else:
                    # Polyfit log(tau) and n
                    # We pad to 6 coefficients (5th degree polynomial)
                    z1 = np.polyfit(T[indices], np.log(tau[indices]), 5)
                    z2 = np.polyfit(T[indices], n[indices], 0)

                    # Ensure consistent 6-length coefficient arrays
                    z1 = np.pad(z1, (6 - len(z1), 0), 'constant')
                    z2 = np.pad(z2, (6 - len(z2), 0), 'constant')
            else:
                # KM Logic (Martensite)
                ms, beta = ttt_data[0:2]
                z1 = np.nan_to_num([ms], nan=273.15)
                z2 = np.nan_to_num([beta], nan=0.01)

            Z1_list.append(z1)
            Z2_list.append(z2)

        # Convert to arrays for interpolation
        X_points = np.array(X_list)
        Z1_points = np.array(Z1_list)
        Z2_points = np.array(Z2_list)

        # 4. Multi-dimensional Interpolation setup
        # Create axis coordinates for interpn
        inter_axes = []
        for i in range(len(comp_elements)):
            axis = sorted(list(set(X_points[:, i])))
            inter_axes.append(axis)
        inter_axes = [sorted(list(set(X_points[:, i]))) for i in range(len(comp_elements))]

        # Reshape Z data to match the meshgrid of inter_axes
        grid_shape = [len(ax) for ax in inter_axes]
        grid_shape = tuple(len(ax) for ax in inter_axes)

        # Prepare the grid of mesh nodes for interpolation
        target_grid = np.stack([full_composition[el] for el in comp_elements], axis=-1)

        print(f"Interpolating {phase} parameters...")

        # Perform interpolation for each polynomial coefficient
        res1_cols, res2_cols = [], []
        num_coeffs = Z1_points.shape[1]
        from scipy.interpolate import RBFInterpolator
        for i in range(num_coeffs):
            # Reshape values into the multi-dim composition grid
            #v1 = Z1_points[:, i].reshape(grid_shape)
            #v2 = Z2_points[:, i].reshape(grid_shape)

            #r1 = interpolate.interpn(inter_axes, v1, target_grid, method="linear", bounds_error=False, fill_value=None)
            #r2 = interpolate.interpn(inter_axes, v2, target_grid, method="linear", bounds_error=False, fill_value=None)
            # 'linear' or 'thin_plate_spline' are usually safe bets
            # neighbors=None uses all points (accurate but slower for huge datasets)
            interp1 = RBFInterpolator(X_points, Z1_points[:, i], kernel='linear')
            interp2 = RBFInterpolator(X_points, Z2_points[:, i], kernel='linear')

            r1 = interp1(target_grid)
            r2 = interp2(target_grid)

            res1_cols.append(r1)
            res2_cols.append(r2)

        # Stack results (N_nodes, N_coeffs)
        final_res1 = np.stack(res1_cols, axis=-1)
        final_res2 = np.stack(res2_cols, axis=-1)

        #import matplotlib.pyplot as plt
        #plt.plot(target_grid[:,0], np.array(final_res1),'o')
        #print(full_composition)
        #plt.figure()
        #plt.plot(full_composition["C"],'o')
        #plt.show()

        # 5. Save to Datastream
        prefix = "JMAK" if phase != "Martensite" else "KM"
        names = ["tau", "n"] if phase != "Martensite" else ["Ms", "b"]

        adjustdatastream({f"{prefix}_{names[0]}_{phase}": final_res1}, datapos="nodes")
        adjustdatastream({f"{prefix}_{names[1]}_{phase}": final_res2}, datapos="nodes")

    print("Success: Transformation models interpolated to all nodes.")


def TTTpolyfit_old():
    from scipy import interpolate
    data = read_geninput()

    # Getting composition for all nodes
    fullcomposition = dict()
    for element in data['Material']['Composition'].keys():
        fullcomposition[element] = readdatastream("Composition_" + element)

    # Getting all compositions from database
    compositions = getTTTcompositions()
    newcomp = list()
    for comp in compositions:
        if comp not in newcomp:
            newcomp.append(comp)
    compositions = newcomp.copy()


    print("Number of compositions used for interpolations are " + str(len(compositions)))

    # Starting polyfit
    for phase in ["Ferrite", "Bainite", "Pearlite", "Martensite"]:
        # Dummy names for JMAK and KM parameters
        Z1 = ""
        Z2 = ""
        #
        X = []
        gridlen = len(data['Material']['Composition'])
        for comp in compositions:
            x = list()
            TTTdata = getTTTdb(comp, "Modeldata")

            if phase in ["Ferrite", "Bainite", "Pearlite"]:
                T = TTTdata[phase][0]
                tau = TTTdata[phase][1]  # tau
                n = TTTdata[phase][2]  # n
                # Checking the values for tau and n. If all nan, n = 1, tau = 1E12
                taumax = 1E8
                if np.isnan(n).all():
                    n = np.nan_to_num(n, nan=3.)
                else:
                    n = np.nan_to_num(n, nan=3.)
                    #n = np.nan_to_num(n, nan=float(np.nanmean(n)))
                tau = np.nan_to_num(tau, nan=taumax)
                tau[tau > taumax] = taumax

                indx = [i for i, v in enumerate(tau) if v < taumax]
                # print(indx)
                polynomial = 4

                if len(indx) == 0:
                    z1 = np.concatenate([np.zeros(5), [taumax]])
                    z2 = np.concatenate([np.zeros(5), [3.]])
                else:
                    indx = [np.min(indx) - 1] + indx
                    if np.max(indx) < len(T)-1:
                        indx = indx + [np.max(indx) + 1]

                    import matplotlib.pyplot as plt
                    # z2 = [0.3]
                    #z2 = [0.5]
                    z1 = np.polyfit(T[indx], np.log(tau[indx]), polynomial)
                    z2 = np.polyfit(T[indx], n[indx], 0)
                    #plt.plot(T[indx], tau[indx])
                    #plt.plot(T[indx], np.exp(np.polyval(z1, T[indx])))
                    #plt.yscale('log')
                    #plt.show()

                    # z2 = [0.4]
                    #z2 = [0.001]
                    # print(z2)
                    while len(z1) < 6:
                        z1 = [0] + list(z1)
                    while len(z2) < 6:
                        z2 = [0] + list(z2)
                    z1 = np.array(z1)
                    z2 = np.array(z2)
                                    # CHANGED TO 4 INSTEAD OF 2
                    # z1 = np.polyfit(T[indx], tau[indx], 2)
                    # z2 = np.polyfit(T[indx], n[indx], 2)

                    # ptest = np.poly1d(z2)
                    # plt.plot(ptest(T[indx]), T[indx])
                    # plt.plot(np.log(n[indx]), T[indx])
                    # plt.show()
            else:  # Martensite
                z1 = [TTTdata[phase][0]]  # Ms
                z2 = [TTTdata[phase][1]]  # beta
                # Checking the values for tau and n if all nan, Ms = 0, beta = 0.01
                z1 = np.nan_to_num(z1, nan=273.15)
                z2 = np.nan_to_num(z2, nan=0.01)

            # Creating new grid with points
            for element in data["Material"]["Composition"].keys():
                tmpx = [comp[element]]
                x.append(tmpx)

            x = np.transpose(x)

            if len(X) != 0:
                X = np.vstack([X, x])
            else:
                X = x.copy()
            # TTT data points
            if isinstance(Z1, str):
                Z1 = z1.copy()
                Z2 = z2.copy()
            elif phase in ["Ferrite", "Bainite", "Pearlite"]:
                Z1 = np.vstack([Z1, z1])
                Z2 = np.vstack([Z2, z2])
            else:
                Z1 = np.vstack([Z1, z1])
                Z2 = np.vstack([Z2, z2])
        # X is all points on the [T, comp1, comp2, comp3, ...] list
        # Z is all points n, tau, Ms, beta values

        # Adjusting grid for interpolation
        tmpX = list(set(np.array(X)[:, 0]))
        tmpX.sort()
        interX = [tmpX]
        for elementnr in range(gridlen-1):
            tmpX = list(set(np.array(X)[:, elementnr+1]))
            tmpX.sort()
            interX.append(tmpX)

        # Interpolating
        print("Interpolating modeldata to gridpoints for " + phase)
        res1 = ""
        res2 = ""
        print(np.shape(Z1))
        try:
            looping = range(len(Z1[0]))
        except:
            looping = [0]
        for i in looping:
            # Adding data to all points and making sure the position is correlating to composition
            interZ1 = np.transpose(np.array(Z1[:, i]).reshape(np.shape(np.transpose(np.meshgrid(*interX, indexing='ij')[0]))))
            interZ2 = np.transpose(np.array(Z2[:, i]).reshape(np.shape(np.transpose(np.meshgrid(*interX, indexing='ij')[0]))))

            # Creating array with composition for all meshpoints
            grid = list()
            for element in data["Material"]["Composition"].keys():
                if len(grid) == 0:
                    grid = [[j] for j in np.array(fullcomposition[element])]
                else:
                    tmpgrid = np.array(fullcomposition[element])
                    grid = [grid[j] + [tmpgrid[j]] for j in range(len(tmpgrid))]


            print("Interpolating polynom parameter " + str(i+1) + "/" + str(len(Z1[0])))

            # Rounding the compositions for easier handling and making sure everything is within interpolation range
            for j in range(len(grid)):
                grid[j] = [round(grid[j][k], 4) for k in range(2)] + [round(grid[j][k], 1) for k in range(2, len(grid[j]))]

            # Interpolating between data points and grid points
            z1 = list(interpolate.interpn(interX, interZ1, grid, method="linear"))
            z2 = list(interpolate.interpn(interX, interZ2, grid, method="linear"))

            # Checking if res is empty or not. If empty copy first result.
            if isinstance(res1, str):
                res1 = z1.copy()
                res2 = z2.copy()
            else:
                res1 = np.vstack([res1, z1])
                res2 = np.vstack([res2, z2])

        # Transposing result before adding to datastream
        if len(res1) == len(Z1[0]):
            res1 = res1.transpose()
            res2 = res2.transpose()


        if phase in ["Ferrite", "Bainite", "Pearlite"]:
            adjustdatastream({"JMAK_tau_" + phase: np.array(res1)})
            adjustdatastream({"JMAK_n_" + phase: np.array(res2)})
        else:
            adjustdatastream({"KM_Ms_" + phase: np.array(res1)})
            adjustdatastream({"KM_b_" + phase: np.array(res2)})
        print(phase + " transformation model added to datastream")
    print("Transformation models interpolated to nodes")
