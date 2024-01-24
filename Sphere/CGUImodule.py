import matplotlib.backends.backend_tkagg
from HelpFile import *
import sys
import threading
import logging
import matplotlib as mpl
from matplotlib.figure import Figure
import customtkinter as ctk
from queue import Queue


class PrintLogger(object):
    def __init__(self, textbox, log_level=logging.INFO):  # pass reference to text widget
        self.textbox = textbox  # keep ref
        self.log_level = log_level

    def write(self, text):
        self.textbox.configure(state="normal")  # make field editable
        self.textbox.insert("end", text)  # write text to textbox
        self.textbox.see("end")  # scroll to end
        self.textbox.configure(state="disabled")  # make field readonly
        self.textbox.update()
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
        from PIL import Image
        #self.columnconfigure(0, weight=1)

        bearingimage = ctk.CTkImage(light_image=Image.open("GUIfiles/Ballbearing2.png"), size=(50, 60))
        imagelabel = ctk.CTkLabel(self, text="", image=bearingimage)
        imagelabel.grid(row=0, column=0, sticky="w", pady=(20, 10), padx=(20, 0))

        self.logo_label = ctk.CTkLabel(self, text="Quenching simulation",
                                       font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=1, padx=20, pady=(20, 10), sticky="w")
        #self.logo_label = ctk.CTkLabel(self, text="Material is: ",
        #                               font=ctk.CTkFont(size=20, weight="bold"))
        #self.logo_label.grid(row=0, column=1, padx=20, pady=(20, 10))


class leftFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.log_widget = ctk.CTkTextbox(self,
                                         width=400,
                                         height=800)
        self.log_widget.grid(row=0, column=0, columnspan=2, padx=20, pady=10, sticky="nsew")

        #self.runall_switch = ctk.CTkSwitch(self, text="Run all modules")
        #self.runall_switch.grid(row=1, column=0, padx=20, pady=20)
        self.sidebar_button_1 = ctk.CTkButton(self, text="Continue")
        self.sidebar_button_1.grid(row=1, column=1, padx=20, pady=20)

        self.progress_bar = ctk.CTkProgressBar(self, height=40,
                                               corner_radius=0,
                                               mode="determinate",
                                               progress_color="green")
        self.progress_bar.set(0)
        self.progress_bar.grid(row=1, column=0, padx=20, pady=20, sticky="nsw")

        #self.sidebar_button_2 = ctk.CTkButton(self, text="Test")
        #self.sidebar_button_2.grid(row=1, column=0, padx=20, pady=20)

    def test(self, master):
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
        tab0frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
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
        if type == "Meshing":
            tab1 = self.tabs_frame.add("Mesh")
            tab1frame = meshTab(tab1)
            tab1frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
            tab1.rowconfigure(0, weight=1)
            tab1.columnconfigure(0, weight=1)
            self.tabs_frame.set("Mesh")
        elif type == "Carbonitriding":
            tab2 = self.tabs_frame.add("Carbonitriding")
            tab2frame = CNTab(tab2)
            tab2frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
            tab2.rowconfigure(0, weight=1)
            tab2.columnconfigure(0, weight=1)
            self.tabs_frame.set("Carbonitriding")
        elif type == "TTT":
            tab3 = self.tabs_frame.add("TTT diagrams")
            tab3frame = TTTTab(tab3)
            tab3frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
            tab3.rowconfigure(0, weight=1)
            tab3.columnconfigure(0, weight=1)
            self.tabs_frame.set("TTT diagrams")
        elif type == "TTTmodeling":
            tab4 = self.tabs_frame.add("TTTmodel")
            tab4frame = TTTmodelTab(tab4)
            tab4frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
            tab4.rowconfigure(0, weight=1)
            tab4.columnconfigure(0, weight=1)
            self.tabs_frame.set("TTTmodel")
        elif type == "Quenching":
            tab5 = self.tabs_frame.add("Quenching")
            tab5frame = QuenchingTab(tab5)
            tab5frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
            tab5.rowconfigure(0, weight=1)
            tab5.columnconfigure(0, weight=1)
            self.tabs_frame.set("Quenching")
        elif type == "End":
            pass


class infoTab(ctk.CTkScrollableFrame):
    def __init__(self, master):
        super().__init__(master, fg_color="transparent")
        import matplotlib as mpl
        from PIL import Image
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        data = read_input()
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)


        # bearingimage = ctk.CTkImage(light_image=Image.open("GUIfiles/Ballbearing.png"),size=(300,300))
        # imagelabel = ctk.CTkLabel(self, text="", image=bearingimage)
        # imagelabel.grid(row=0, column=0, sticky="nsew")
        compstr = [key +  ": " + str(data["Material"]["Composition"][key]) for key in data["Material"]["Composition"].keys()]
        infostrings = ["Composition is:    " + "    ".join(compstr),
                       "Number of nodes in radius: " + str(data["Geometry"]["nodes"]),
                       "Radius of sphere: " + str(data["Geometry"]["radius"]),
                       "Mesh program: " + data["Programs"]["Meshing"],
                       "Carbonitriding program: " + data["Programs"]["Carbonitriding"],
                       "TTT program: " + data["Programs"]["TTT"],
                       "FEM program: " + data["Programs"]["FEM"]]
        i = 0
        for info in infostrings:
            self.rowconfigure(i, weight=0)
            self.composition = ctk.CTkLabel(self, text=info, pady=10)
            self.composition.grid(row=i, column=0, sticky="w")
            i = i + 1

        mpl.rcParams["font.size"] = 32
        text = ctk.CTkLabel(self, text="Temperature history:")
        text.grid(row=i, column=0, columnspan=2, sticky="nsew")

        i = i + 1
        # Add plot of temperature
        starttemp = data["Thermo"]["CNtemp"] - 273.15
        quenchtemp = data["Thermo"]["quenchtemp"] - 273.15
        tempertemp = 400  # data["Thermo"]["tempertemp"]
        roomtemp = data["Thermo"]["quenchtemp"] - 273.15
        holdCN = data["Thermo"]["CNtime"]
        holdquench = holdCN + 1800
        holdtemper = holdquench + 1800
        holdend = holdtemper + 1800
        times = np.array([0, 0, holdCN, holdCN + 1800, holdquench, holdquench + 1800, holdtemper, holdtemper + 1800, holdend]) / 3600
        temps = [roomtemp, starttemp, starttemp, quenchtemp, quenchtemp, tempertemp, tempertemp, roomtemp, roomtemp]
        fig, ax = plt.subplots(figsize=(6, 10), dpi=50)
        ax.plot(times, temps)
        ax.set_xlabel('Time [h]')
        ax.set_ylabel('Temperature [degC]')

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew", padx=20)
        # toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        # toolbar.grid(row=5, column=0, sticky="nsew")
        # toolbar.update()
        canvas._tkcanvas.grid(row=i, column=0, columnspan=2, sticky="nsew")
        self.rowconfigure(i, weight=1)


class meshTab(ctk.CTkFrame):
    def __init__(self, master):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        super().__init__(master, fg_color="transparent")
        mpl.rcParams["font.size"] = 32
        self.columnconfigure(0, weight=1)
        # self.rowconfigure(0, weight=1)
        # self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        grid = meshio.read("Resultfiles/Datastream.xdmf")
        nodes = grid.points

        self.nodenr = ctk.CTkLabel(self, text="Number of nodes: " + str(len(nodes)), pady=10)
        self.nodenr.grid(row=0, column=0, sticky="nsew")
        # self.choices = ctk.CTkComboBox(self, values=["test 1", "test 2"])
        # self.choices.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")
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
    def __init__(self, master, fg_color="transparent"):
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
        plot1.set_xlabel('Radius [mm]')
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
    def __init__(self, master):
        super().__init__(master, fg_color="transparent")
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
        super().__init__(master, fg_color="transparent")
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

class FEMresults(Figure):
    def __init__(self, master, datatype):
        super().__init__((10, 4), 50, fg_color="transparent")
        self.master = master
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        cord = getaxisvalues('nodes')
        if datatype == "vonMises":
            azm = np.linspace(0, 2 * np.pi, 30)
            r, th = np.meshgrid(cord[:, 0], azm)
            vonMises = [getaxisvalues("vonMises") for i in range(30)]
            ax = self.add_subplot(projection="polar")
            plot3 = ax.pcolormesh(th, r, vonMises)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xlabel("testing")
            self.colorbar(plot3, ax=self.get_axes(), label="MPa")

            self.canvas = FigureCanvasTkAgg(self, master=master)
            self.canvas.draw()
            self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
            toolbar = NavigationToolbar2Tk(self.canvas, master, pack_toolbar=False)
            toolbar.grid(row=2, column=0, sticky="nsew")
            toolbar.update()
            self.tkcanvas.grid(row=1, column=0, sticky="nsew")

    def show(self):
        pass

    def clear(self):
        pass


class QuenchingTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        cord = getaxisvalues('nodes')

        def combobox_callback(choice):
            widgets = self.grid_slaves(row=1, column=0)
            if len(widgets) != 0:
                for widget in widgets:
                    widget.grid_remove()
            canvas[choice]._tkcanvas.grid(row=1, column=0, sticky="nsew")

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

        fig.colorbar(plot1, ax=subfigs[0].get_axes(), label="Austenite fraction")
        fig.colorbar(plot2, ax=subfigs[1].get_axes(), label="Martensite fraction")
        plot1.set_clim(0, 1)
        plot2.set_clim(0, 1)

        phasecanvas = FigureCanvasTkAgg(fig, master=self)
        phasecanvas.draw()
        phasecanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(phasecanvas, self, pack_toolbar=False)
        toolbar.grid(row=2, column=0, sticky="nsew")
        toolbar.update()
        phasecanvas._tkcanvas.grid(row=0, column=0, sticky="nsew")

        fig2 = Figure(figsize=(10, 4), dpi=50)
        vonMises = [getaxisvalues("vonMises") for i in range(30)]

        ax3 = fig2.add_subplot(projection="polar")
        plot3 = ax3.pcolormesh(th, r, vonMises)
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])

        fig2.colorbar(plot3, ax=fig2.get_axes())

        stresscanvas = FigureCanvasTkAgg(fig2, master=self)
        stresscanvas.draw()
        stresscanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(stresscanvas, self, pack_toolbar=False)
        toolbar.grid(row=2, column=0, sticky="nsew")
        toolbar.update()
        # stresscanvas._tkcanvas.grid(row=1, column=0, sticky="nsew")

        strainfigure = Figure(figsize=(10, 4), dpi=50)
        ep1 = getaxisvalues("ep1")
        ep2 = getaxisvalues("ep2")
        ep3 = getaxisvalues("ep3")
        strainaxis = strainfigure.add_subplot()
        strainaxis.plot(cord[:, 0], ep1, label="First principal strain")
        strainaxis.plot(cord[:, 0], ep2, label="Second principal strain")
        strainaxis.plot(cord[:, 0], ep3, label="Third principal strain")
        strainaxis.legend()
        strainaxis.set_xticklabels([])
        strainaxis.set_yticklabels([])
        strainaxis.set_xlabel('Radius [mm]')
        strainaxis.set_ylabel('Strain [-]')

        straincanvas = FigureCanvasTkAgg(strainfigure, master=self)
        straincanvas.draw()
        straincanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(straincanvas, self, pack_toolbar=False)
        toolbar.grid(row=2, column=0, sticky="nsew")
        toolbar.update()

        canvas = dict()
        canvas["von Mises stress"] = stresscanvas
        canvas["Phase fractions"] = phasecanvas
        canvas["Principal strain"] = straincanvas
        # test = FEMresults(self, "vonMises")
        #test.show()
        combobox = ctk.CTkComboBox(master=self,
                                   values=list(canvas.keys()),
                                   command=combobox_callback)
        combobox.grid(row=0, column=0, sticky="nsew")
        # print("combobox dropdown clicked:", choice)


class MainApp(ctk.CTk):
    def __init__(self):
        from StructureFile import CalcModule
        super().__init__()
        self.geometry("1400x1000")
        self.title("Quenching of steel")

        # Setting weights
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        # Create header
        self.header_frame = headerFrame(self)
        self.header_frame.grid(row=0, column=0, sticky="nsew")

        # Create sidebar
        self.sidebar_frame = leftFrame(self)
        self.sidebar_frame.grid(row=1, column=0, columnspan=1, sticky="nsew")

        # Creating resultwindow
        self.main_frame = rightFrame(self)
        self.main_frame.grid(row=0, column=1, rowspan=2, columnspan=1, sticky="nsew")

        # Adding button functionality
        #self.sidebar_frame.sidebar_button_1.configure(command=lambda: threading.Thread(target=self.next_module).start())
        self.sidebar_frame.sidebar_button_1.configure(command=self.next_module)
        #self.sidebar_frame.sidebar_button_2.configure(command=lambda: threading.Thread(target=self.test).start())

        self.programstate = ctk.IntVar(self, 0)
        #self.runall = ctk.IntVar(self, self.sidebar_frame.runall_switch.get())
        self.modules = Queue()
        self.modules.put(CalcModule("Meshing"))
        self.modules.put(CalcModule("Carbonitriding"))
        self.modules.put(CalcModule("TTT"))
        self.modules.put(CalcModule("TTTmodeling"))
        #self.modules.put(CalcModule("Quenching"))
    def test(self):
        logger = PrintLogger(self.sidebar_frame.log_widget)
        sys.stdout = logger
        sys.stderr = logger
        import time
        print("Testing module")
        for i in range(1, 3):
            time.sleep(5)
            print(str(i * 5) + "sec")

    def next_module(self):
        self.sidebar_frame.progress_bar.set(0.0)
        logger = PrintLogger(self.sidebar_frame.log_widget)
        sys.stdout = logger
        #sys.stderr = logger

        if self.programstate.get() == 0:
            createdatastreamcache()
            print("Created a cache file from previous results\n")

        if self.modules.empty():
            print("\nNo modules left in pipeline")
            self.sidebar_frame.sidebar_button_1.grid_remove()
            self.sidebar_frame.progress_bar.grid_remove()
            savedatastream("Resultfiles/230124_1.xdmf")
            print("Resultdata saved to " + "Resultfiles/230124_1.xdmf")
            return

        currentmodule = self.modules.get()
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=self.run_module, args=(currentmodule,))
            #tid.daemon = True
            tid.start()
            self.progressmonitor(tid, currentmodule)

        else:
            tid = None
            self.run_module(currentmodule)
            self.sidebar_frame.progress_bar.set(currentmodule.getprogress())
        self.programstate.set(self.programstate.get() + 1)
        data = read_input()
    def run_module(self, module):
        module.runmodule()
        self.main_frame.add_gui(module.modulename())

    def progressmonitor(self, tid, module):
        #print("Progress is " + str(module.getprogress()))
        self.sidebar_frame.progress_bar.set(module.getprogress())
        """ Monitor the download thread """
        if tid.is_alive():
            self.after(100, lambda: self.sidebar_frame.progress_bar.set(module.getprogress()))
            pass
            #self.after(1000, lambda: self.progressmonitor(tid, module))
            #self.after(100, lambda: self.sidebar_frame.progress_bar.set(module.getprogress()))
