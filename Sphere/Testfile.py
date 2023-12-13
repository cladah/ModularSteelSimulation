import sqlite3

import numpy as np
from sqlitedict import SqliteDict
from HelpFile import *
import matplotlib.pyplot as plt
#db_file = 'Resultfiles/database.db'
#with sqlite3.connect(db_file) as conn:
#    print('Created the connection!')
#composition = SqliteDict("Resultfiles/database.db", tablename="composition", outer_stack=False)

def test1():
    import h5py
    import numpy as np

    analyseTTTdatabase()
    input("pause")



    data = read_input()
    with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
        x = np.array(f.get("CNcurves/Position"))
        a = dict()
        for element in data['Material']['Composition'].keys():
            #plt.plot(x, 100 * np.array(f.get("CNcurves/"+element)))
            a[element] = round(100 * np.array(f.get("CNcurves/Elements/" + element))[-1], 2) # Composition curves at all points along x

    #print(a)

    filename = "TTT_center.hdf5"

    phases = ["Ferrite","Bainite","Perlite","Martensite"]
    TTTdata = dict()
    modeldata = dict()
    for ph in phases:
        phase = dict()
        model = dict()
        if ph in ["Ferrite", "Bainite", "Perlite"]:
            Tsteps = readresultfile(filename, ph + "/Tsteps")
            phase["start"] = [Tsteps, readresultfile(filename, ph + "/start")]
            phase["half"] = [Tsteps, readresultfile(filename, ph + "/half")]
            phase["finish"] = [Tsteps, readresultfile(filename, ph + "/finish")]
            model["JMAK"] = [readresultfile(filename,ph + "/JMAK/T"),
                             readresultfile(filename,ph +"/JMAK/tau"),
                             readresultfile(filename,ph +"/JMAK/n")]
        else:
            start = readresultfile(filename, ph + "/start")
            half = readresultfile(filename, ph + "/half")
            finish = readresultfile(filename, ph + "/finish")
            phase["start"] = [[start,start],[0.1, 1E12]]
            phase["half"] = [[half,half],[0.1, 1E12]]
            phase["finish"] = [[finish,finish],[0.1, 1E12]]
            model["KM"] = [readresultfile(filename, ph + "/KM/Ms"),
                           readresultfile(filename, ph + "/KM/beta")]

        TTTdata[ph] = phase
        modeldata[ph] = model
    data = read_input()
    composition = data["Material"]["Composition"]
    #addTTTdata(composition, TTTdata, "TTTdata")
    #addTTTdata(composition, modeldata, "Modeldata")
    #e = getTTTdata(composition,"TTTdata")
    #print(e)
    #TTTdata = SqliteDict("Resultfiles/database.db", tablename="TTTdata", outer_stack=False)
    #TTTdata["1"] = {"TTTdata":1,"Composition":2}
    #print([i for i in TTTdata.keys()])
    #for key, item in TTTdata.items():
    #    print(key)
    #print(TTTdata.items())
    #TTTdata.close()
    #test = getTTTdata(composition)
    #plt.plot(test["Perlite"]["start"][1],test["Perlite"]["start"][0])
    #plt.plot(test["Bainite"]["start"][1],test["Bainite"]["start"][0])
    #plt.plot(test["Ferrite"]["start"][1],test["Ferrite"]["start"][0])
    #plt.plot(test["Martensite"]["start"][1],test["Martensite"]["start"][0])
    #plt.plot(test["Martensite"]["finish"][1],test["Martensite"]["finish"][0])
    #plt.xscale("log")
    #plt.axis([1,1E12,300,1000])
    #plt.show()


    #print(readdatastream("C"))

def splineJMAK(phase, filename):
    import numpy as np
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

def modeldatatocsv():
    import csv
    xyz = readdatastream("nodes")
    phases = ["Ferrite", "Perlite", "Bainite"]
    for phase in phases:
        if phase != "Martensite":
            tmpT = readresultfile("Modeldata", phase + "/JMAK/T")
            tmp1 = readresultfile("Modeldata", phase + "/JMAK/tau")
            tmp2 = readresultfile("Modeldata", phase + "/JMAK/n")
            tmpnames = ["tau", "n"]
        else:
            tmpT = readresultfile("Modeldata", phase + "/T")
            tmp1 = readresultfile("Modeldata", phase + "/Ms")
            tmp2 = readresultfile("Modeldata", phase + "/beta")
            tmpnames = ["Ms", "beta"]

        for i in [0, 1]:
        #     with open("Resultfiles/" + phase + "_" + tmpnames[i] + ".txt", 'w') as f:
        #         if not tmp1.all():
        #             break
        #         if i == 0:
        #             for j in range(len(tmp1)-1):
        #                 for k in range(len(tmpT)-1):
        #                     f.write(" ".join(np.str([xyz[j][0],xyz[j][1],tmpT[k],tmp1[j][k]])) + "\n")
        #         else:
        #             for j in range(len(tmp2)-1):
        #                 for k in range(len(tmpT) - 1):
        #                     f.write(" ".join([xyz[j][0],xyz[j][1],tmpT[k],tmp2[j][k]]) + "\n")
            with open("Resultfiles/" + phase + "_" + tmpnames[i] + ".csv", 'w') as file:
             ##with open("Resultfiles/" + phase + "_" + tmpnames[i] + ".txt", 'w') as file:
                 writer = csv.writer(file)
                 if not tmp1.all():
                     break
                 if i == 0:
                     for j in range(len(tmp1)-1):
                         for k in range(len(tmpT)-1):
                             writer.writerow([xyz[j][0], xyz[j][1], tmpT[k], tmp1[j][k]])
                 else:
                     for j in range(len(tmp2)-1):
                         for k in range(len(tmpT) - 1):
                             writer.writerow([xyz[j][0], xyz[j][1], tmpT[k], tmp2[j][k]])
        print(str(phase) + " written to csv files")