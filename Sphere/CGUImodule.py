from HelpFile import *
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitridingmodule
from TTTmodule import runTTTmodule, runTTTmodelmodule
from Quenchingmodule import runquenchingmodule
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
class testFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.grid_columnconfigure(0, weight=1)
        self.logo_label = ctk.CTkLabel(self, text="Quenching simulation",
                                       font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(column=0, padx=20, pady=(20, 10))
class headerFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.columnconfigure(1, weight=1)

        self.logo_label = ctk.CTkLabel(self, text="Quenching simulation",
                                                 font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.logo_label = ctk.CTkLabel(self, text="Material is: ",
                                       font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=1, padx=20, pady=(20, 10))
class leftFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)


        self.log_widget = ctk.CTkTextbox(self,
                                         width=400,
                                         height=800)
        self.log_widget.grid(row=0, column=0, columnspan=2, padx=20, pady=10, sticky="nsew")

        self.sidebar_button_1 = ctk.CTkButton(self, text="Continue")
        self.sidebar_button_1.grid(row=1, column=1, padx=20, pady=20)
    def test(self,master):
        print()
class rightFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        # Adding a tab
        self.tabs_frame = ctk.CTkTabview(self)
        self.tabs_frame.grid(row=1, column=0, padx=20, pady=10, sticky="nsew")
        tab0 = self.tabs_frame.add("Input")
        tab0frame = infoTab(tab0)
        tab0frame.grid(row=0, column=0, sticky="nsew")
        tab0.rowconfigure(0, weight=1)
        tab0.columnconfigure(0, weight=1)
    def add_gui(self, type):
        from HelpFile import readdatastream, getaxisvalues, getTTTdata, read_input
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        mpl.rcParams["font.size"] = 32
        if type == "Mesh":
            tab1 = self.tabs_frame.add("Mesh")
            tab1frame = meshTab(tab1)
            tab1frame.grid(row=0, column=0, sticky="nsew")
            tab1.rowconfigure(0, weight=1)
            tab1.columnconfigure(0, weight=1)
            self.tabs_frame.set("Mesh")
        elif type == "Carbonitriding":
            tab2 = self.tabs_frame.add("Carbonitriding")
            tab2frame = CNTab(tab2)
            tab2frame.grid(row=0, column=0, sticky="nsew")
            tab2.rowconfigure(0, weight=1)
            tab2.columnconfigure(0, weight=1)
            self.tabs_frame.set("Carbonitriding")
        elif type == "TTT":
            tab3 = self.tabs_frame.add("TTT diagrams")
            tab3frame = TTTTab(tab3)
            tab3frame.grid(row=0, column=0, sticky="nsew")
            tab3.rowconfigure(0, weight=1)
            tab3.columnconfigure(0, weight=1)
            self.tabs_frame.set("TTT diagrams")
        elif type == "TTTmodel":
            tab4 = self.tabs_frame.add("TTTmodel")
            tab4frame = TTTmodelTab(tab4)
            tab4frame.grid(row=0, column=0, sticky="nsew")
            tab4.rowconfigure(0, weight=1)
            tab4.columnconfigure(0, weight=1)
            self.tabs_frame.set("TTTmodel")
        elif type == "Quenching":
            tab5 = self.tabs_frame.add("Quenching")
            tab5frame = QuenchingTab(tab5)
            tab5frame.grid(row=0, column=0, sticky="nsew")
            tab5.rowconfigure(0, weight=1)
            tab5.columnconfigure(0, weight=1)
            self.tabs_frame.set("Quenching")
        elif type == "End":
            pass
class infoTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        data = read_input()
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)
        self.nodenr = ctk.CTkLabel(self, text="Number of nodes: ", pady=10)
        self.nodenr.grid(row=0, column=0, sticky="nsew")

        data = read_input()
        mpl.rcParams["font.size"] = 32

        text = ctk.CTkLabel(self, text="Temperature history:")
        text.grid(row=1, column=0, sticky="nsew")

        # Add plot of temperature
        starttemp = data["Thermo"]["CNtemp"] - 273.15
        quenchtemp = data["Thermo"]["quenchtemp"] - 273.15
        tempertemp = 400  # data["Thermo"]["tempertemp"]
        roomtemp = data["Thermo"]["quenchtemp"] - 273.15
        holdCN = data["Thermo"]["CNtime"] * 60
        holdquench = holdCN + 30
        holdtemper = holdquench + 60
        holdend = holdtemper + 60
        times = [0, 0, holdCN, holdCN + 10, holdquench, holdquench + 10, holdtemper, holdtemper + 10, holdend]
        temps = [roomtemp, starttemp, starttemp, quenchtemp, quenchtemp, tempertemp, tempertemp, roomtemp, roomtemp]
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        plot1.plot(times, temps)
        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=2, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=3, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=2, column=0, sticky="nsew")
class meshTab(ctk.CTkFrame):
    def __init__(self, master):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        super().__init__(master)
        mpl.rcParams["font.size"] = 32
        self.columnconfigure(0, weight=1)
        #self.rowconfigure(0, weight=1)
        #self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        grid = meshio.read("Resultfiles/Datastream.xdmf")
        nodes = grid.points

        self.nodenr = ctk.CTkLabel(self, text="Number of nodes: " + str(len(nodes)), pady=10)
        self.nodenr.grid(row=0, column=0, sticky="nsew")
        #self.choices = ctk.CTkComboBox(self, values=["test 1", "test 2"])
        #self.choices.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")
        self.elnr = ctk.CTkLabel(self, text="Number of elements: " + str(len(nodes)), pady=0)
        self.elnr.grid(row=1, column=0, sticky="nsew")

        fig, ax = plt.subplots(figsize=(5, 4), dpi=50)
        cmap = mpl.colors.ListedColormap("lightgray")
        c = np.ones(len(nodes))
        ax.tripcolor(nodes[:, 0] * 1000, nodes[:, 1] * 1000, c, edgecolor="k", cmap=cmap)
        fig.gca().set_aspect('equal')
        ax.set_xlabel('[mm]')
        ax.set_ylabel('[mm]')

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=2, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=3, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=2, column=0, sticky="nsew")
class CNTab(ctk.CTkFrame):
    def __init__(self,master):
        super().__init__(master)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

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
        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")

class TTTTab(ctk.CTkFrame):
    def __init__(self,master):
        super().__init__(master)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

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

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")

class TTTmodelTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        cord = getaxisvalues('nodes')
        fig = Figure(figsize=(10, 4), dpi=50)

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")
class QuenchingTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
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

        phasecanvas = FigureCanvasTkAgg(fig, master=self)
        phasecanvas.draw()
        phasecanvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(phasecanvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        #phasecanvas._tkcanvas.grid(row=0, column=0, sticky="nsew")

        fig2 = Figure(figsize=(10, 4), dpi=50)
        vonMises = [getaxisvalues("vonMises") for i in range(30)]

        ax3 = fig2.add_subplot(projection="polar")
        plot3 = ax3.pcolormesh(th, r, vonMises)
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])

        fig2.colorbar(plot3, ax=fig2.get_axes())


        stresscanvas = FigureCanvasTkAgg(fig2, master=self)
        stresscanvas.draw()
        stresscanvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(stresscanvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        stresscanvas._tkcanvas.grid(row=1, column=0, sticky="nsew")

        canvas = dict()
        canvas["Stress"] = stresscanvas
        canvas["Phase fraction"] = phasecanvas

        def combobox_callback(choice):
            self.grid_slaves(row=1,column=0)[0].grid_remove()
            canvas[choice]._tkcanvas.grid(row=1, column=0, sticky="nsew")
            # print("combobox dropdown clicked:", choice)

        combobox = ctk.CTkComboBox(master=self,
                                             values=["Phase fractions", "von Mises stress"],
                                             command=combobox_callback)

        combobox.grid(row=0, column=0, sticky="nsew")


class MainApp(ctk.CTk):
    def __init__(self):
        from StructureFile import CalcModule
        super().__init__()
        self.geometry("1200x1000")
        self.title("Quenching of steel")

        # Setting weights
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        # Create header
        self.header_frame = headerFrame(self)
        self.header_frame.grid(row=0, column=0, columnspan=2, sticky="nsew")

        # Create sidebar
        self.sidebar_frame = leftFrame(self)
        self.sidebar_frame.grid(row=1, column=0, columnspan=1, sticky="nsew")


        # Creating resultwindow
        self.main_frame = rightFrame(self)
        self.main_frame.grid(row=1, column=1, columnspan=1, sticky="nsew")

        # Adding button functionality
        self.sidebar_frame.sidebar_button_1.configure(command=self.next_module)

        # self.sidebar_frame = ctk.CTkFrame(self, width=300, corner_radius=0)
        # self.sidebar_frame.grid(row=0, column=3, rowspan=4, sticky="nsew")
        # self.sidebar_frame.grid_rowconfigure(4, weight=1)
        # self.sidebar_button_1 = ctk.CTkButton(self.sidebar_frame, command=self.next_module)
        # self.sidebar_button_1.grid(row=1, column=0, padx=20, pady=10)

        # terminal
        #self.log_widget = ctk.CTkScrollableFrame.ScrolledText(self.sidebar_frame, height=50, width=50,
        #                                               font=("consolas", "12", "normal"))
        #self.log_widget.pack()


        # Create tabframe
        #self.main_frame = ctk.CTkFrame(self, corner_radius=0)
        #self.main_frame.grid(row=0, column=1, sticky="nsew")
        #self.tabs_frame = ctk.CTkTabview(self.main_frame)
        #self.tabs_frame.grid(row=0, column=0, sticky="nsew")
        #self.sidebar_button_1 = ctk.CTkButton(self.tabs_frame, command=self.next_module)
        #self.sidebar_button_1.grid(row=1, column=0, padx=20, pady=10)

        self.programstate = ctk.IntVar(self, 0)
        self.modules = list()
        self.modules.append(CalcModule("Meshing"))
        self.modules.append(CalcModule("Carbonitriding"))
        self.modules.append(CalcModule("TTT"))
        self.modules.append(CalcModule("TTTmodeling"))
        self.modules.append(CalcModule("Quenching"))

    def next_module(self):
        # Logging frame
        logger = PrintLogger(self.sidebar_frame.log_widget)
        sys.stdout = logger
        sys.stderr = logger
        modules = ["Mesh","Carbonitriding","TTT","TTTmodel","Quenching","Tempering"]

        # GUI frame
        try:
            threading.Thread(self.run_module(self.programstate.get())).start()
            self.main_frame.add_gui(modules[self.programstate.get()])
        except IndexError:
            print("No modules left in pipeline")
        #self.tabs.select(self.programstate.get())
        self.programstate.set(self.programstate.get() + 1)

    def run_module(self, type):
        self.modules[type].runmodule()
        # if type == 0:
        #
        # elif type == 1:
        #     runcarbonitridingmodule()
        # elif type == 2:
        #     runTTTmodule()
        # elif type == 3:
        #     runTTTmodelmodule()
        # elif type == 4:
        #     runquenchingmodule()