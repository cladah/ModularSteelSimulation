from HelpFile import *
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitridingmodule
from TTTmodule import runTTTmodule, runTTTmodelmodule
from Quenchingmodule import runquenchingmodule
import tkinter as tk
from tkinter import scrolledtext
from tkinter import ttk
import sys
import threading
from multiprocessing import Process, Pipe
import subprocess
from io import StringIO
import logging
import matplotlib as mpl
import customtkinter as ctk


class PrintLogger(object):
    def __init__(self, textbox, log_level=logging.INFO):  # pass reference to text widget
        self.textbox = textbox  # keep ref
        self.log_level = log_level

    def write(self, text):
        self.textbox.configure(state="normal")  # make field editable
        self.textbox.insert("end", text)  # write text to textbox
        self.textbox.see("end")  # scroll to end
        self.textbox.configure(state="disabled")  # make field readonly
    def flush(self):  # needed for file like object
        pass




class MainApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.geometry("1200x1000")
        self.title("Quenching of steel")
        self.root = tk.Frame(self)
        self.root.pack()
        self.programstate = tk.IntVar(self, 0)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(1, weight=1)

        self.header = tk.Frame(self.root, width=100, height=200, bg="blue")
        self.header.grid(row=0, columnspan=2, sticky="news")
        self.header.columnconfigure(0, weight=1)
        self.header.rowconfigure(0, weight=1)



        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        self.left_frame = tk.Frame(self.root, width=200, height=200, bg='grey')
        self.left_frame.grid(column=0, row=1, sticky="w", columnspan=1)
        self.left_frame.columnconfigure(0, weight=1)
        self.left_frame.rowconfigure(1, weight=1)

        self.right_frame = tk.Frame(self.root, width=200, height=300, bg='grey')
        self.right_frame.grid(column=1, row=1, sticky="e", columnspan=1)
        self.right_frame.columnconfigure(1, weight=1)
        self.right_frame.rowconfigure(1, weight=1)
        button = tk.Button(master=self.left_frame, text="Continue", command=self.next_module, padx=20, pady=10)
        button.pack(side=tk.TOP, anchor=tk.NW, padx=20, pady=20)
        self.log_widget = tk.scrolledtext.ScrolledText(self.left_frame, height=50, width=50,
                                                       font=("consolas", "12", "normal"))
        self.log_widget.pack()
        #self.left_frame.pack(expand=False)
        #self.right_frame.pack(expand=False)
        #textbox = tk.Text(self.right_frame, height=10, width=80, wrap='word')
        #textbox.pack(side=tk.TOP, anchor=tk.NW, padx=20, pady=20)
        self.tabs = ttk.Notebook(self.right_frame)
        self.tabs.pack(side="top", fill="both", expand=True)
        self.start_page()


    def start_page(self):
        data = read_input()
        mpl.rcParams["font.size"] = 32
        tab = ttk.Frame(self)

        text = tk.Label(tab, text="Temperature history:")
        text.pack()
        self.tabs.pack()
    def reset_logging(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    def next_module(self):
        logger = PrintLogger(self.log_widget)
        sys.stdout = logger
        sys.stderr = logger
        modules = ["Mesh","Carbonitriding","TTT","TTTmodel","Quenching","Tempering"]


        threading.Thread(self.run_module(self.programstate.get())).start()
        self.add_gui(modules[self.programstate.get()])
        self.tabs.select(self.programstate.get())
        self.programstate.set(self.programstate.get() + 1)

    def run_module(self, type):
        if type == 0:
            createMesh()
        elif type == 1:
            runcarbonitridingmodule()
        elif type == 2:
            runTTTmodule()
        elif type == 3:
            runTTTmodelmodule()
        elif type == 4:
            runquenchingmodule()

    def add_gui(self, type):
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

            tab = ttk.Frame(self)

            nodenr = tk.Label(tab, text="Number of nodes: " + str(len(nodes)), pady=10)
            nodenr.pack()
            elnr = tk.Label(tab, text="Number of elements: " + str(len(nodes)), pady=0)
            elnr.pack()

            fig, ax = plt.subplots(figsize=(5, 4), dpi=50)
            cmap = mpl.colors.ListedColormap("lightgray")
            c = np.ones(len(nodes))
            ax.tripcolor(nodes[:, 0] * 1000, nodes[:, 1] * 1000, c, edgecolor="k", cmap=cmap)
            fig.gca().set_aspect('equal')
            ax.set_xlabel('[mm]')
            ax.set_ylabel('[mm]')

            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            canvas.get_tk_widget().pack()
            toolbar = NavigationToolbar2Tk(canvas)
            toolbar.update()
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            self.tabs.add(tab, text="Mesh")
            self.tabs.pack()

        elif type == "Carbonitriding":
            tab = ttk.Frame(self)
            wC = getaxisvalues("Composition/C")
            wN = getaxisvalues("Composition/N")
            xyz = getaxisvalues("nodes")
            fig = Figure(figsize=(5, 4), dpi=50)
            plot1 = fig.add_subplot(111)
            plot1.plot(np.array(xyz)[:, 0] * 1000, wC, label='Carbon')
            plot1.plot(np.array(xyz)[:, 0] * 1000, wN, label='Nitrogen')
            plot1.set_xlabel('radius [mm]')
            plot1.set_ylabel('Weight fraction [%]')
            plot1.legend()
            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            canvas.get_tk_widget().pack()
            toolbar = NavigationToolbar2Tk(canvas)
            toolbar.update()
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            self.tabs.add(tab, text="Carbonitriding")
            self.tabs.pack()

        elif type == "TTT":
            tab = ttk.Frame(self)
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
            colorlist = ["green", "blue", "orange", "red"]
            i = 0
            for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
                plot1.plot(TTTcore[phase]["start"][1], TTTcore[phase]["start"][0], label=phase, color=colorlist[i])
                plot1.plot(TTTcore[phase]["finish"][1], TTTcore[phase]["finish"][0], linestyle="dashed",
                           color=colorlist[i])
                i = i + 1
            plot1.set_xscale('log')
            plot1.title.set_text('Core TTT')
            plot1.legend(loc="upper right")
            plot2 = fig.add_subplot(122)
            plot2.set_xlim([0.1, 1.E12])
            i = 0
            for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
                plot2.plot(TTTsurf[phase]["start"][1], TTTsurf[phase]["start"][0], label=phase, color=colorlist[i])
                plot2.plot(TTTsurf[phase]["finish"][1], TTTsurf[phase]["finish"][0], linestyle="dashed",
                           color=colorlist[i])
                i = i + 1
            plot2.set_xscale('log')
            plot2.title.set_text('Surface TTT')
            plot2.legend(loc="upper right")
            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            canvas.get_tk_widget().pack()
            self.tabs.add(tab, text="TTT diagrams")
            self.tabs.pack()
        elif type == "TTTmodel":
            tab = ttk.Frame(self)
            self.tabs.add(tab, text="TTT fitted")
            self.tabs.pack()
        elif type == "Quenching":
            tab = ttk.Frame(self)
            cord = getaxisvalues('nodes')
            fig = Figure(figsize=(10, 4), dpi=50)

            azm = np.linspace(0, 2 * np.pi, 30)
            r, th = np.meshgrid(cord[:, 0], azm)
            aust = [getaxisvalues("Austenite") for i in range(30)]
            mart = [getaxisvalues("Martensite") for i in range(30)]

            # sfig1 = fig.add_subplot(121, projection="polar")
            # sfig2 = fig.add_subplot(122, projection="polar")

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
            self.tabs.add(tab, text="Results")
            self.tabs.pack()

    def initialize_logging(self):
        logger = PrintLogger(self.log_widget)
        sys.stdout = logger
        sys.stderr = logger

        print("Program has started")