import matplotlib.pyplot
import numpy as np

from HelpFile import *
import sys
import threading
import logging
from matplotlib.figure import Figure
import customtkinter as ctk
from queue import Queue
import meshio
from Datastream_file import getaxisvalues, readdatastream, createdatastreamcache, savedatastream, getnamesdatastream, gethistoryvalues
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule
from Postprocessing.Postprocess_main import read_input_result

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
        import matplotlib as mpl
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
        elif type == "Transformationmodels":
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
            #self.tabs_frame.tab("TTTmodel")
        elif type == "End":
            pass


class infoTab(ctk.CTkScrollableFrame):
    def __init__(self, master):
        super().__init__(master, fg_color="transparent")
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
        data = read_input()
        self.columnconfigure(0, weight=0)
        self.columnconfigure(3, weight=1)


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
            self.composition.grid(row=i, column=0, columnspan=4, sticky="w")
            i = i + 1
        self.rowconfigure(i, weight=0)
        self.composition = ctk.CTkLabel(self, text="Cachefile", pady=10)
        self.composition.grid(row=i, column=0, columnspan=1, sticky="w")
        self.refbox = ctk.CTkTextbox(self, pady=2, width=400, height=30, wrap="none")
        self.refbox.grid(row=i, column=1, columnspan=1, sticky="nsew", padx=10)
        self.refbox.insert("end", text=data["Datastream"]["Cachedirect"])


        def savesettings():
            logger = PrintLogger(self.master.master.master.master.master.master.sidebar_frame.log_widget)
            logger.write("Added " + self.refbox.get("1.0",ctk.END).strip("\n") + " as cache file.\n")
            change_input("Datastream", "Cachedirect", self.refbox.get("1.0",ctk.END).strip("\n"))
            self.refbox.configure(state=ctk.DISABLED)
            self.save_button.configure(state=ctk.DISABLED)
            self.master.master.master.master.master.master.input = read_input_result(self.refbox.get("1.0",ctk.END).strip("\n"))
        self.save_button = ctk.CTkButton(self, text="Load", command=savesettings)
        self.save_button.grid(row=i, column=2, columnspan=1, sticky="nsew", padx=10)





        i = i + 1

        # mpl.rcParams["font.size"] = 32
        # text = ctk.CTkLabel(self, text="Temperature history:")
        # text.grid(row=i, column=0, columnspan=4, sticky="nsew")
        #
        # i = i + 1
        # Add plot of temperature
        starttemp = data["Thermo"]["CNtemp"] - 273.15
        quenchtemp = data["Thermo"]["quenchtemp"] - 273.15
        tempertemp = 400  # data["Thermo"]["tempertemp"]
        roomtemp = data["Thermo"]["quenchtemp"] - 273.15
        holdCN = data["Thermo"]["CNtime"]
        holdquench = holdCN + 1800
        holdtemper = holdquench + 1800
        holdend = holdtemper + 1800
        plt.rcParams.update({'font.size': 22})
        times = np.array([0, 0, holdCN, holdCN + 600, holdquench, holdquench + 1800, holdend]) / 3600
        temps = [roomtemp, starttemp, starttemp, quenchtemp, quenchtemp, roomtemp, roomtemp]
        fig, ax = plt.subplots(figsize=(6, 10), dpi=50)
        ax.plot(times, temps)
        ax.set_xlabel('Time [h]')
        ax.set_ylabel('Temperature [degC]')

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew", padx=20, pady=20)
        # toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        # toolbar.grid(row=5, column=0, sticky="nsew")
        # toolbar.update()
        canvas._tkcanvas.grid(row=i, column=0, columnspan=4, sticky="nsew")
        self.rowconfigure(i, weight=1)


class meshTab(ctk.CTkFrame):
    def __init__(self, master):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        super().__init__(master, fg_color="transparent")
        mpl.rcParams["font.size"] = 32
        self.columnconfigure(0, weight=1)
        # self.rowconfigure(0, weight=1)
        # self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        nodes = readdatastream("nodes")
        elements = readdatastream("elements")

        # print(elements)
        # print(elements[0])
        # print(len(elements[0]))

        self.nodenr = ctk.CTkLabel(self, text="Number of nodes: " + str(len(nodes)), pady=10)
        self.nodenr.grid(row=0, column=0, sticky="nsew")
        # self.choices = ctk.CTkComboBox(self, values=["test 1", "test 2"])
        # self.choices.grid(row=0, column=0, padx=20, pady=10, sticky="nsew")
        self.elnr = ctk.CTkLabel(self, text="Number of elements: " + str(len(elements[0])), pady=0)
        self.elnr.grid(row=1, column=0, sticky="nsew")

        fig, ax = plt.subplots(figsize=(5, 4), dpi=50)

        yz = np.c_[nodes[:, 0] * 1000, nodes[:, 1]* 1000]
        #print(elements[0].data)

        cmap = mpl.colors.ListedColormap("lightgray")
        # reordering vert
        elemdata = elements[0].data
        elemdata = [e[[0, 3, 1, 4, 2, 5]] for e in elemdata]
        #print(elemdata)
        #verts = yz[elements[0].data]
        verts = yz[elemdata]
        #print(verts)
        pc = mpl.collections.PolyCollection(verts, color="lightgray", edgecolor="k")
        #print(pc)
        ax.add_collection(pc)


        c = np.ones(len(nodes))
        #ax.tripcolor(nodes[:, 0] * 1000, nodes[:, 1] * 1000, c, edgecolor="k", cmap=cmap)
        fig.gca().set_aspect('equal')
        ax.set_xlabel('[mm]')
        ax.set_ylabel('[mm]')
        ax.set_xlim([-np.max(nodes[:, 0]*1000)*0.1, np.max(nodes[:, 0]*1000)*1.1])
        ax.set_ylim([-np.max(nodes[:, 1]*1000)*0.1, np.max(nodes[:, 1]*1000)*1.1])


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
        plot1.set_ylabel('Weight [%]')
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

        fig = Figure(figsize=(20, 8), dpi=50)
        plot1 = fig.add_subplot(121)
        plot1.set_xlim([0.1, 1.E12])
        colorlist = ["green", "blue", "orange", "red"]
        i = 0
        for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
            plot1.plot(TTTcore[phase]["start"][1], np.array(TTTcore[phase]["start"][0])-273.15, label=phase, color=colorlist[i])
            plot1.plot(TTTcore[phase]["finish"][1], np.array(TTTcore[phase]["finish"][0])-273.15, linestyle="dashed",
                       color=colorlist[i])
            i = i + 1
        plot1.grid()
        plot1.set_xscale('log')
        plot1.title.set_text('Core TTT')
        plot1.legend(loc="upper right")
        plot1.set_xlabel('Time [s]')
        plot1.set_ylabel('Temperature [degC]')
        plot1.set_ylim([0, 900])

        plot2 = fig.add_subplot(122)
        plot2.set_xlim([0.1, 1.E12])
        i = 0
        for phase in ["Ferrite", "Bainite", "Perlite", "Martensite"]:
            plot2.plot(TTTsurf[phase]["start"][1], np.array(TTTsurf[phase]["start"][0])-273.15, label=phase, color=colorlist[i])
            plot2.plot(TTTsurf[phase]["finish"][1], np.array(TTTsurf[phase]["finish"][0])-273.15, linestyle="dashed",
                       color=colorlist[i])
            i = i + 1
        plot2.grid()
        plot2.set_xscale('log')
        plot2.title.set_text('Surface TTT')
        plot2.set_xlabel('Time [s]')
        plot2.set_ylabel('Temperature [degC]')
        plot2.legend(loc="upper right")
        plot2.set_ylim([0, 900])
        fig.tight_layout()
        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        self.canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")


class TTTmodelTab(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master, fg_color="transparent")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        Tgrid = np.linspace(0,1000,100)
        fig = Figure(figsize=(20, 8), dpi=50)
        plot1 = fig.add_subplot(121)
        plot1.set_xlim([0.1, 1.E12])
        #colorlist = ["green", "blue", "orange", "red"]
        colorlist = ["blue", "orange", "red"]
        i = 0
        for phase in ["Bainite", "Perlite", "Martensite"]:
            if phase in ["Ferrite", "Bainite", "Perlite"]:
                z1 = getaxisvalues("JMAK_tau_" + phase)[0]
                z2 = getaxisvalues("JMAK_n_" + phase)[0]
                p1 = np.poly1d(z1)
                p2 = np.poly1d(z2)
                # Z98 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.02)) ** np.array(np.exp(p2(Tgrid)))
                # Z02 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.98)) ** np.array(np.exp(p2(Tgrid)))
                Z98 = np.array(np.exp(p1(Tgrid))*10) * (-np.log(0.02)) ** (1/np.array(p2(Tgrid)))
                Z02 = np.array(np.exp(p1(Tgrid))*10) * (-np.log(0.98)) ** (1/np.array(p2(Tgrid)))
                indx = [i for i, v in enumerate(Z98) if v < 1E12]
                Z98 = Z98[indx]
                Z02 = Z02[indx]
                X = Tgrid[indx]
                plot1.plot(Z02, X - 273.15, label=phase,
                           color=colorlist[i])
                plot1.plot(Z98, X - 273.15, linestyle="dashed",
                           color=colorlist[i])

            else:
                z1 = getaxisvalues("KM_Ms_" + phase)[0]
                z2 = getaxisvalues("KM_b_" + phase)[0]
                start = z1 + np.log(0.98) / z2 - 273.15
                finish = z1 + np.log(0.02) / z2 - 273.15
                plot1.plot([0.1, 1E12], [start, start], label=phase,
                           color=colorlist[i])
                plot1.plot([0.1, 1E12], [finish, finish], linestyle="dashed",
                           color=colorlist[i])
            i = i + 1
        plot1.grid()
        plot1.set_xscale('log')
        plot1.title.set_text('Core TTT')
        plot1.legend(loc="upper right")
        plot1.set_xlabel('Time [s]')
        plot1.set_ylabel('Temperature [degC]')
        plot1.set_ylim([0, 900])

        plot2 = fig.add_subplot(122)
        plot2.set_xlim([0.1, 1.E12])
        i = 0
        for phase in ["Bainite", "Perlite", "Martensite"]:
            if phase in ["Ferrite", "Bainite", "Perlite"]:
                z1 = getaxisvalues("JMAK_tau_" + phase)[-1]
                z2 = getaxisvalues("JMAK_n_" + phase)[-1]
                p1 = np.poly1d(z1)
                p2 = np.poly1d(z2)
                Z98 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.02)) ** (1/np.array(p2(Tgrid)))
                Z02 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.98)) ** (1/np.array(p2(Tgrid)))
                indx = [i for i, v in enumerate(Z98) if v < 1E12]
                Z98 = Z98[indx]
                Z02 = Z02[indx]
                X = Tgrid[indx]
                plot2.plot(Z02, X - 273.15, label=phase,
                           color=colorlist[i])
                plot2.plot(Z98, X - 273.15, linestyle="dashed",
                           color=colorlist[i])

            else:
                z1 = getaxisvalues("KM_Ms_" + phase)[-1]
                z2 = getaxisvalues("KM_b_" + phase)[-1]
                start = z1 + np.log(0.98)/z2 - 273.15
                finish = z1 + np.log(0.02)/z2 - 273.15
                plot2.plot([0.1, 1E12], [start, start], label=phase,
                           color=colorlist[i])
                plot2.plot([0.1,1E12], [finish, finish], linestyle="dashed",
                           color=colorlist[i])
            i = i + 1
        plot2.grid()
        plot2.set_xscale('log')
        plot2.title.set_text('Surface TTT')
        plot2.set_xlabel('Time [s]')
        plot2.set_ylabel('Temperature [degC]')
        plot2.legend(loc="upper right")
        plot2.set_ylim([0, 900])
        fig.tight_layout()
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
        #self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)
        #self.rowconfigure(0, weight=1)

        def combobox_callback(choice):
            widgets = self.grid_slaves(row=1, column=0)
            if len(widgets) != 0:
                for widget in widgets:
                    widget.grid_remove()
            canvas[choice]._tkcanvas.grid(row=1, column=0, sticky="nsew")
            #canvas[choice].get_tk_widget().grid(row=1, column=0, sticky="nsew")
            toolbar = NavigationToolbar2Tk(canvas[choice], self, pack_toolbar=False)
            toolbar.grid(row=2, column=0, sticky="nsew")
            toolbar.update()

        aust = getaxisvalues("Austenite", time=-1)
        mart = getaxisvalues("Martensite", time=-1)
        fer = getaxisvalues("Ferrite", time=-1)
        per = getaxisvalues("Perlite", time=-1)
        bai = getaxisvalues("Bainite", time=-1)
        xyz = getaxisvalues("nodes")
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        plot1.plot(np.array(xyz)[:, 0] * 1000, aust, label="Austenite", color="purple")
        #plot1.plot(np.array(xyz)[:, 0] * 1000, fer, label="Ferrite", color="green")
        plot1.plot(np.array(xyz)[:, 0] * 1000, per, label="Perlite", color="orange")
        plot1.plot(np.array(xyz)[:, 0] * 1000, bai, label="Bainite", color="blue")
        plot1.plot(np.array(xyz)[:, 0] * 1000, mart, label="Martensite", color="red")
        plot1.set_xlabel('Radius [mm]')
        plot1.set_ylabel('Phase fraction [-]')
        plot1.legend()
        plot1.title.set_text('Phase fraction after 600s')
        phasecanvas = FigureCanvasTkAgg(fig, master=self)
        phasecanvas.draw()
        phasecanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(phasecanvas, self, pack_toolbar=False)
        toolbar.grid(row=2, column=0, sticky="nsew")
        toolbar.update()
        phasecanvas._tkcanvas.grid(row=1, column=0, sticky="nsew")

        s1 = getaxisvalues("vonMises", time=-1)
        s = getaxisvalues("Stress", time=-1)
        vM2 = np.sqrt(s[:,0]**2+s[:,2]**2+s[:,5]**2 - s[:,2]*s[:,5]-s[:,0]*s[:,2]-s[:,0]*s[:,5] - 3*(s[:,1]**2 + s[:,3]**2+s[:,4]**2))
        #sp = getaxisvalues("P_Stress", time=-1)
        #sp1 = sp[:, 0]
        #sp2 = sp[:, 1]
        #sp3 = sp[:, 2]
        #R = ((stress[:,0]-stress[:,2])/2)**2+
        # s2 = getaxisvalues("sl22", time=-1)
        # s3 = getaxisvalues("sl33", time=-1)
        xyz = getaxisvalues("nodes")
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        #plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(s1)*1e-6, label="von-Mises stress")
        plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(vM2) * 1e-6, label="von-Mises stress")
        # plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(sp1) * 1e-6, label="First principal")
        # plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(sp2) * 1e-6, label="Second principal")
        # plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(sp3) * 1e-6, label="Third principal")
        # plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(s2)*1e-6, label="Second principal")
        # plot1.plot(np.array(xyz)[:, 0] * 1000, np.array(s3)*1e-6, label="Third principal")
        plot1.set_xlabel('Radius [mm]')
        plot1.set_ylabel('Stress [MPa]')
        plot1.legend()
        plot1.title.set_text('vonMises stress after 600s')
        stresscanvas = FigureCanvasTkAgg(fig, master=self)
        stresscanvas.draw()
        stresscanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(phasecanvas, self, pack_toolbar=False)
        toolbar.grid(row=2, column=0, sticky="nsew")
        toolbar.update()
        #stresscanvas._tkcanvas.grid(row=0, column=0, sticky="nsew")


        # TTT with temp
        Tgrid = np.linspace(0, 1000, 100)
        fig = Figure(figsize=(20, 8), dpi=50)
        plot1 = fig.add_subplot(121)
        plot1.set_xlim([0.1, 1.E12])
        colorlist = ["blue", "orange", "red"]
        i = 0
        for phase in ["Bainite", "Perlite", "Martensite"]:
            if phase in ["Ferrite", "Bainite", "Perlite"]:
                z1 = getaxisvalues("JMAK_tau_" + phase)[0]
                z2 = getaxisvalues("JMAK_n_" + phase)[0]
                p1 = np.poly1d(z1)
                p2 = np.poly1d(z2)
                Z98 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.02)) ** np.array(p2(Tgrid))
                Z02 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.98)) ** np.array(p2(Tgrid))
                indx = [j for j, v in enumerate(Z98) if v < 1E12]
                Z98 = Z98[indx]
                Z02 = Z02[indx]
                X = Tgrid[indx]
                plot1.plot(Z02, X - 273.15,
                           color=colorlist[i])
                plot1.plot(Z98, X - 273.15, linestyle="dashed",
                           color=colorlist[i])
            else:
                z1 = getaxisvalues("KM_Ms_" + phase)[0]
                z2 = getaxisvalues("KM_b_" + phase)[0]
                start = z1 + np.log(0.98)/z2 - 273.15
                finish = z1 + np.log(0.02)/z2 - 273.15
                plot1.plot([0.1, 1E12], [start, start],
                           color=colorlist[i])
                plot1.plot([0.1,1E12], [finish, finish], linestyle="dashed",
                           color=colorlist[i])
            i = i + 1
        timex, tempcore = gethistoryvalues("T", 0)
        plot1.plot(timex, np.array(tempcore) - 273.15, label="Temperature core",
                   color="black")
        plot1.grid()
        plot1.set_xscale('log')
        plot1.title.set_text('Core TTT')
        plot1.set_xlabel('Time [s]')
        plot1.set_ylabel('Temperature [degC]')
        plot1.legend(loc="upper right")
        plot1.set_ylim([0, 900])

        plot2 = fig.add_subplot(122)
        plot2.set_xlim([0.1, 1.E12])
        colorlist = ["blue", "orange", "red"]
        i = 0
        for phase in ["Bainite", "Perlite", "Martensite"]:
            if phase in ["Ferrite", "Bainite", "Perlite"]:
                z1 = getaxisvalues("JMAK_tau_" + phase)[-1]
                z2 = getaxisvalues("JMAK_n_" + phase)[-1]
                p1 = np.poly1d(z1)
                p2 = np.poly1d(z2)
                Z98 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.02)) ** np.array(p2(Tgrid))
                Z02 = np.array(np.exp(p1(Tgrid))) * (-np.log(0.98)) ** np.array(p2(Tgrid))
                indx = [j for j, v in enumerate(Z98) if v < 1E12]
                Z98 = Z98[indx]
                Z02 = Z02[indx]
                X = Tgrid[indx]
                plot2.plot(Z02, X - 273.15,
                           color=colorlist[i])
                plot2.plot(Z98, X - 273.15, linestyle="dashed",
                           color=colorlist[i])
            else:
                z1 = getaxisvalues("KM_Ms_" + phase)[-1]
                z2 = getaxisvalues("KM_b_" + phase)[-1]
                start = z1 + np.log(0.98)/z2 - 273.15
                finish = z1 + np.log(0.02)/z2 - 273.15
                plot2.plot([0.1, 1E12], [start, start],
                           color=colorlist[i])
                plot2.plot([0.1,1E12], [finish, finish], linestyle="dashed",
                           color=colorlist[i])
            i = i + 1
        timex, tempsurf = gethistoryvalues("T", -1)
        plot2.plot(timex, np.array(tempsurf) - 273.15, label="Temperature surface", linestyle="dashed",
                   color="black")
        plot2.grid()
        plot2.set_xscale('log')
        plot2.title.set_text('Surface TTT')
        plot2.set_xlabel('Time [s]')
        plot2.set_ylabel('Temperature [degC]')
        plot2.legend(loc="upper right")
        plot2.set_ylim([0, 900])

        fig.tight_layout()
        TTTcanvas = FigureCanvasTkAgg(fig, master=self)
        TTTcanvas.draw()
        TTTcanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(TTTcanvas, self, pack_toolbar=False)
        toolbar.grid(row=2, column=0, sticky="nsew")
        toolbar.update()


        # def combobox_callback(choice):
        #     widgets = self.grid_slaves(row=1, column=0)
        #     if len(widgets) != 0:
        #         for widget in widgets:
        #             widget.grid_remove()
        #     canvas[choice]._tkcanvas.grid(row=1, column=0, sticky="nsew")
        #
        # fig = Figure(figsize=(10, 4), dpi=50)
        #
        # azm = np.linspace(0, 2 * np.pi, 30)
        # r, th = np.meshgrid(cord[:, 0], azm)
        #
        # aust = [getaxisvalues("Austenite", time=-1) for i in range(30)]
        # mart = [getaxisvalues("Martensite", time=-1) for i in range(30)]
        #
        # # sfig1 = fig.add_subplot(121, projection="polar")
        # # sfig2 = fig.add_subplot(122, projection="polar")
        #
        # subfigs = fig.subfigures(1, 2)
        #
        # ax1 = subfigs[0].add_subplot(projection="polar")
        # plot1 = ax1.pcolormesh(th, r, aust)
        # ax1.set_xticklabels([])
        # ax1.set_yticklabels([])
        #
        # ax2 = subfigs[1].add_subplot(projection="polar")
        # plot2 = ax2.pcolormesh(th, r, mart)
        # ax2.set_xticklabels([])
        # ax2.set_yticklabels([])
        #
        # fig.colorbar(plot1, ax=subfigs[0].get_axes(), label="Austenite fraction")
        # fig.colorbar(plot2, ax=subfigs[1].get_axes(), label="Martensite fraction")
        # plot1.set_clim(0, 1)
        # plot2.set_clim(0, 1)
        #
        # phasecanvas = FigureCanvasTkAgg(fig, master=self)
        # phasecanvas.draw()
        # phasecanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        # toolbar = NavigationToolbar2Tk(phasecanvas, self, pack_toolbar=False)
        # toolbar.grid(row=2, column=0, sticky="nsew")
        # toolbar.update()
        # phasecanvas._tkcanvas.grid(row=0, column=0, sticky="nsew")
        #
        # fig2 = Figure(figsize=(10, 4), dpi=50)
        # vonMises = [getaxisvalues("vonMises") for i in range(30)]
        #
        # ax3 = fig2.add_subplot(projection="polar")
        # plot3 = ax3.pcolormesh(th, r, vonMises)
        # ax3.set_xticklabels([])
        # ax3.set_yticklabels([])
        #
        # fig2.colorbar(plot3, ax=fig2.get_axes())
        #
        # stresscanvas = FigureCanvasTkAgg(fig2, master=self)
        # stresscanvas.draw()
        # stresscanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        # toolbar = NavigationToolbar2Tk(stresscanvas, self, pack_toolbar=False)
        # toolbar.grid(row=2, column=0, sticky="nsew")
        # toolbar.update()
        # # stresscanvas._tkcanvas.grid(row=1, column=0, sticky="nsew")
        #
        # strainfigure = Figure(figsize=(10, 4), dpi=50)
        # ep1 = getaxisvalues("ep1",time=-1)
        # ep2 = getaxisvalues("ep2",time=-1)
        # ep3 = getaxisvalues("ep3",time=-1)
        # strainaxis = strainfigure.add_subplot()
        # strainaxis.plot(cord[:, 0], ep1, label="First principal strain")
        # strainaxis.plot(cord[:, 0], ep2, label="Second principal strain")
        # strainaxis.plot(cord[:, 0], ep3, label="Third principal strain")
        # strainaxis.legend()
        # strainaxis.set_xticklabels([])
        # strainaxis.set_yticklabels([])
        # strainaxis.set_xlabel('Radius [mm]')
        # strainaxis.set_ylabel('Strain [-]')
        #
        # straincanvas = FigureCanvasTkAgg(strainfigure, master=self)
        # straincanvas.draw()
        # straincanvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        # toolbar = NavigationToolbar2Tk(straincanvas, self, pack_toolbar=False)
        # toolbar.grid(row=2, column=0, sticky="nsew")
        # toolbar.update()
        #
        canvas = dict()
        canvas["Phase fractions"] = phasecanvas
        # canvas["Strain"] = straincanvas
        canvas["Stress"] = stresscanvas
        canvas["Temprature history"] = TTTcanvas
        # test = FEMresults(self, "vonMises")
        #test.show()
        combobox = ctk.CTkComboBox(master=self,
                                   values=list(canvas.keys()),
                                   command=combobox_callback)
        combobox.grid(row=0, column=0, pady=(10, 10), sticky="nsew")
        # # print("combobox dropdown clicked:", choice)


class MainApp(ctk.CTk):
    def __init__(self):

        super().__init__()
        self.wm_iconbitmap("GUIfiles/Ballbearing.ico")
        self.geometry("1600x1000")
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
        self.modules.put(Meshingmodule())
        self.modules.put(Carbonitridingmodule())
        self.modules.put(TTTdiagrammodule())
        self.modules.put(Transformationmodelmodule())
        self.modules.put(Quenchingmodule())
        self.input = read_input()


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
        self.sidebar_frame.sidebar_button_1.configure(state=ctk.DISABLED)
        self.sidebar_frame.progress_bar.set(0.0)
        logger = PrintLogger(self.sidebar_frame.log_widget)
        sys.stdout = logger
        #sys.stderr = logger


        if self.programstate.get() == 0:
            self.input = read_input()
            createdatastreamcache(self.input["Datastream"]["Cachedirect"])

        if self.modules.empty():
            print("\nNo modules left in pipeline")
            self.sidebar_frame.sidebar_button_1.grid_remove()
            self.sidebar_frame.progress_bar.grid_remove()
            savedatastream(self.input["Datastream"]["Savedirect"])
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
            self.sidebar_frame.sidebar_button_1.configure(state=ctk.NORMAL)
        self.programstate.set(self.programstate.get() + 1)
        data = read_input()
    def run_module(self, module):
        module.run()
        self.main_frame.add_gui(module.modulename())

    def progressmonitor(self, tid, module):
        #print("Progress is " + str(module.getprogress()))
        self.sidebar_frame.progress_bar.set(module.getprogress())
        """ Monitor the download thread """
        if tid.is_alive():
            self.after(1000, lambda: self.progressmonitor(tid, module))
        else:
            self.sidebar_frame.sidebar_button_1.configure(state=ctk.NORMAL)
