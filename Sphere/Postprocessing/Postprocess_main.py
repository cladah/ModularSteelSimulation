import meshio
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from Sphere.HelpFile import read_input
import json

def read_input_result(filename):

    f = open('Cachefiles/filename', 'r')
    data = json.load(f)
    f.close()
    return data
def readresultfile(filename, dataname, time=0):
    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            if dataname == "nodes":
                return points
            elif dataname == "elements":
                return cells
            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                if t == time:
                    return point_data[dataname]
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
            if time == -1:
                maxindx = t_list.index(np.max(t_list))
                return pd_list[maxindx][dataname]
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in" + filename + "file. Data that exist is "
                       + str(pd_list[0].keys()))
def getaxisvalues_result(filename, dataname, time=0):
    node_y = readresultfile(filename, 'nodes')[:, 1]
    indx = np.where(node_y == 0)
    data = readresultfile(filename, dataname, time)[indx]
    x = readresultfile(filename, 'nodes')[:, 0][indx]

    # Sorting data in order from x=min(x) to x=max(x)
    indx = np.argsort(x)
    data = np.array(data)[indx]
    return data

def gethistoryvalues_result(filename, dataname, position):
    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()
        if dataname == "nodes":
            return points
        elif dataname == "elements":
            return cells
        pd_list, t_list = list(), list()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            if dataname in point_data.keys():
                t_list.append(t)
                pd_list.append(point_data[dataname][position])
    return t_list, pd_list

def plotcompare(filenamelist, dataname, time=0):
    xvalues, yvalues = list(), list()
    for fn in filenamelist:
        xvalues.append(getaxisvalues_result(fn, 'nodes', time)[:, 0])
        yvalues.append(getaxisvalues_result(fn, dataname, time))
    for i in range(len(xvalues)):
        plt.plot(xvalues[i], yvalues[i],label=str(i))
    plt.legend()
    plt.show()

def plotTTT(filename, quenching=False, position=0):
    plt.rcParams.update({'font.size': 12})
    Tgrid = np.linspace(0, 1000, 100)
    #fig = mpl.figure.Figure(figsize=(10, 4), dpi=50)
    fig, plot1 = plt.subplots()
    #plot1 = fig.add_subplot(111)
    plot1.set_xlim([0.1, 1.E12])
    colorlist = ["green", "blue", "orange", "red"]


    i = 0
    for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
        if phase in ["Ferrite", "Bainite", "Perlite"]:
            z1 = getaxisvalues_result(filename, "JMAK_tau_" + phase)[position]
            z2 = getaxisvalues_result(filename, "JMAK_n_" + phase)[position]
            print(np.shape(getaxisvalues_result(filename, "JMAK_tau_" + phase)))
            p1 = np.poly1d(z1)
            p2 = np.poly1d(z2)
            Z98 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.02)) ** np.array(p2(Tgrid))
            Z02 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.98)) ** np.array(p2(Tgrid))
            indx = [j for j, v in enumerate(Z98) if v < 1E12]
            Z98 = Z98[indx]
            Z02 = Z02[indx]
            X = Tgrid[indx]
            plot1.plot(Z02, X - 273.15, label=phase,
                       color=colorlist[i])
            plot1.plot(Z98, X - 273.15, linestyle="dashed",
                       color=colorlist[i])
        else:
            z1 = getaxisvalues_result(filename, "KM_Ms_" + phase)[position]
            z2 = getaxisvalues_result(filename, "KM_b_" + phase)[position]
            start = z1 + np.log(0.98) / z2 - 273.15
            finish = z1 + np.log(0.02) / z2 - 273.15
            plot1.plot([0.1, 1E12], [start, start], label=phase,
                       color=colorlist[i])
            plot1.plot([0.1, 1E12], [finish, finish], linestyle="dashed",
                       color=colorlist[i])
        i = i + 1
    if quenching:
        timex, tempsurf = gethistoryvalues_result(filename, "T", -1)
        timex, tempcore = gethistoryvalues_result(filename, "T", 0)
        plot1.plot(timex, np.array(tempcore) - 273.15, label="Temperature core",
                   color="black")
        plot1.plot(timex, np.array(tempsurf) - 273.15, label="Temperature surface", linestyle="dashed",
                   color="black")
    plot1.set_xscale('log')
    #plot1.title.set_text('Surface TTT')
    plot1.set_xlabel('Time [s]')
    plot1.set_ylabel('Temperature [degC]')
    plot1.legend(loc="upper right")
    plot1.set_ylim([0, 900])
    # plot1.show()
    #fig.show()
    plt.show()

def plot_stressstrain(filename):
    data = read_input_result(filename)
    fig = Figure(figsize=(10, 4), dpi=50)
    plot1 = fig.add_subplot(111)
    colorlist = ["purple", "green", "blue", "orange", "red"]
    i = 0
    for phase in ["Austenite", "Ferrite", "Perlite", "Bainite", "Martensite"]:
        plot1.plot([0,],[0, ], label=phase, color=colorlist[i])
        i = i + 1
