import tkinter as tk
from tkinter import ttk
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitridingmodule
from TTTmodule import runTTTmodule, runTTTmodelmodule
from HelpFile import *

def runguimodule(gui):
    import tkinter as tk
    #from tkinter import ttk
    from HelpFile import read_input
    import matplotlib as mpl
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                   NavigationToolbar2Tk)
    import threading
    import matplotlib.pyplot as plt
    data = read_input()
    mpl.rcParams["font.size"] = 32
    tab = ttk.Frame(gui)



    text = tk.Label(tab, text="Temperature history:")
    text.pack()

    # Add plot of temperature
    starttemp = data["Thermo"]["CNtemp"]-273.15
    quenchtemp = data["Thermo"]["quenchtemp"]-273.15
    tempertemp = 400 # data["Thermo"]["tempertemp"]
    roomtemp = data["Thermo"]["quenchtemp"]-273.15
    holdCN = data["Thermo"]["CNtime"]*60
    holdquench = holdCN + 30
    holdtemper = holdquench + 60
    holdend = holdtemper + 60
    times = [0, 0, holdCN, holdCN+10, holdquench,holdquench + 10, holdtemper, holdtemper + 10, holdend]
    temps = [roomtemp, starttemp, starttemp, quenchtemp, quenchtemp, tempertemp, tempertemp, roomtemp, roomtemp]
    fig = Figure(figsize=(5, 4), dpi=50)
    plot1 = fig.add_subplot(111)
    plot1.plot(times,temps)
    canvas = FigureCanvasTkAgg(fig, master=tab)
    canvas.draw()
    canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    canvas.get_tk_widget().pack()
    gui.add(tab, text="Heat treatment history")
    gui.pack()

    #return gui

def addgui(gui, type):
    from HelpFile import readdatastream, getaxisvalues, getTTTdata, read_input
    import numpy as np
    import tkinter as tk
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                   NavigationToolbar2Tk)

    mpl.rcParams["font.size"] = 32

    if type == "Mesh":
        grid = meshio.read("Resultfiles/Datastream.xdmf")
        nodes = grid.points

        tab = ttk.Frame(gui)

        nodenr = tk.Label(tab, text="Number of nodes: " + str(len(nodes)), pady=10)
        nodenr.pack()
        elnr = tk.Label(tab, text="Number of elements: " + str(len(nodes)), pady=0)
        elnr.pack()

        fig, ax = plt.subplots(figsize=(5, 4), dpi=50)
        cmap = mpl.colors.ListedColormap("lightgray")
        c = np.ones(len(nodes))
        ax.tripcolor(nodes[:, 0]*1000, nodes[:, 1]*1000, c, edgecolor="k", cmap=cmap)
        fig.gca().set_aspect('equal')
        ax.set_xlabel('[mm]')
        ax.set_ylabel('[mm]')

        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


        gui.add(tab, text="Mesh")
        gui.pack()

    elif type == "Carbonitriding":
        tab = ttk.Frame(gui)
        wC = getaxisvalues("Composition/C")
        wN = getaxisvalues("Composition/N")
        xyz = getaxisvalues("nodes")
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        plot1.plot(np.array(xyz)[:, 0]*1000, wC, label='Carbon')
        plot1.plot(np.array(xyz)[:, 0]*1000, wN, label='Nitrogen')
        plot1.set_xlabel('radius [mm]')
        plot1.set_ylabel('Weight fraction [%]')
        plot1.legend()
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        gui.add(tab, text="Carbonitriding")
        gui.pack()

    elif type == "TTT":
        tab = ttk.Frame(gui)
        data = read_input()
        core = data['Material']['Composition']
        surface = dict()
        for element in data['Material']['Composition'].keys():
            surface[element] = getaxisvalues("Composition/" + element)[-1]
        TTTcore = getTTTdata(core, "TTTdata")
        TTTsurf = getTTTdata(surface, "TTTdata")
        fig = Figure(figsize=(10, 4), dpi=50)
        plot1 = fig.add_subplot(121)
        plot1.set_xlim([0.1, 1.E12])
        colorlist = ["green","blue","orange","red"]
        i = 0
        for phase in ["Ferrite", "Bainite", "Perlite","Martensite"]:
            plot1.plot(TTTcore[phase]["start"][1],TTTcore[phase]["start"][0], label =phase, color=colorlist[i])
            plot1.plot(TTTcore[phase]["finish"][1],TTTcore[phase]["finish"][0], linestyle="dashed", color=colorlist[i])
            i = i + 1
        plot1.set_xscale('log')
        plot1.title.set_text('Core TTT')
        plot1.legend(loc="upper right")
        plot2 = fig.add_subplot(122)
        plot2.set_xlim([0.1, 1.E12])
        i = 0
        for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
            plot2.plot(TTTsurf[phase]["start"][1], TTTsurf[phase]["start"][0], label =phase, color=colorlist[i])
            plot2.plot(TTTsurf[phase]["finish"][1], TTTsurf[phase]["finish"][0], linestyle="dashed", color=colorlist[i])
            i = i + 1
        plot2.set_xscale('log')
        plot2.title.set_text('Surface TTT')
        plot2.legend(loc="upper right")
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        canvas.get_tk_widget().pack()
        gui.add(tab, text="TTT diagrams")
        gui.pack()
    elif type == "TTTmodel":
        tab = ttk.Frame(gui)
        gui.add(tab, text="TTT fitted")
        gui.pack()
    elif type == "Quenching":
        tab = ttk.Frame(gui)
        cord = getaxisvalues('nodes')
        fig = Figure(figsize=(10, 4), dpi=50)

        azm = np.linspace(0, 2 * np.pi, 30)
        r, th = np.meshgrid(cord[:, 0], azm)
        aust = [getaxisvalues("Austenite") for i in range(30)]
        mart = [getaxisvalues("Martensite") for i in range(30)]

        #sfig1 = fig.add_subplot(121, projection="polar")
        #sfig2 = fig.add_subplot(122, projection="polar")

        subfigs = fig.subfigures(1, 2)

        ax1 = subfigs[0].add_subplot(projection="polar")
        plot1 = ax1.pcolormesh(th, r, aust)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])

        ax2 = subfigs[1].add_subplot(projection="polar")
        plot2 = ax2.pcolormesh(th, r, mart)
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])

        fig.colorbar(plot1, ax=subfigs[0].get_axes())
        fig.colorbar(plot2, ax=subfigs[1].get_axes())
        plot1.set_clim(0, 1)
        plot2.set_clim(0, 1)


        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas.draw()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        canvas.get_tk_widget().pack()
        gui.add(tab, text="Results")
        gui.pack()
def _run():
    pass
