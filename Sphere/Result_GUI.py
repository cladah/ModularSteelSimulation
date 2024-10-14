import customtkinter as ctk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
import numpy as np
import matplotlib as mpl
import colorsys
from ResultReading import read_results_history, read_results, getnames_results, read_results_all, read_results_axis
from customtkinter import filedialog

class resultTab(ctk.CTkFrame):
    """
    dataname
    data
    legend
    """


    def __init__(self, master, dataname, data, legdata, xlbl="Time [s]", ylbl = "?"):
        super().__init__(master)
        y = ["Phase fraction [-]", "Phase fraction [-]", "Phase fraction [-]", "Stress [Pa]", "vonMises stress [Pa]",
             "Plastic strain [-]", "Elastic strain [-]", "Temperature [K]",
             "Weight fraction [-]"]
        datanames = ["Austenite", "Bainite", "Martensite", "Stress", "vonMises", "Strain_pl", "Strain", "T"]
        ynames = ["Phase fraction [-]", "Phase fraction [-]", "Phase fraction [-]", "Stress [Pa]", "von-Mises stress [Pa]",
             "Plastic strain [-]", "Elastic strain [-]", "Temperature [K]"]
        strainnames = ["Strain_tot"]

        if "Composition" in dataname:
            ylbl = "Weight [w%]"
            xlbl = "Radius [mm]"
            leg = []
        elif dataname in datanames:
            ylbl = ynames[datanames.index(dataname)]
            xlbl = "Time [s]"
            leg = legdata
        elif dataname == "Phasecomp":
            ylbl = "Phase fraction [-]"
            xlbl = "Radius [mm]"
            leg = legdata
        elif dataname == "Matcomp":
            ylbl = "Weight [w%]"
            xlbl = "Radius [mm]"
            leg = legdata
        elif "Stress" in dataname:
            ylbl = "Stress [Pa]"
            xlbl = "Time [s]"
            leg = legdata
        elif "Strain" in dataname:
            ylbl = "Strain - "
            xlbl = "Time [s]"
            leg = legdata
        else:
            ylbl = "?"
            xlbl = "Time [s]"
            leg = legdata

        def colorFader(c1, c2, mix=0):  # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
            from PIL import ImageColor
            c1 = np.array(ImageColor.getcolor(c1, "RGB"))/255
            c2 = np.array(ImageColor.getcolor(c2, "RGB"))/255
            c1 = np.array(colorsys.rgb_to_hsv(*c1))
            c2 = np.array(colorsys.rgb_to_hsv(*c2))
            c3 = (1. - mix) * c1 + mix * c2
            c3 = tuple(round(i * 255) for i in colorsys.hsv_to_rgb(*c3))
            return "#" + ('{:02X}' * 3).format(*c3)

        mixspace = np.linspace(0., 1., np.shape(data[1])[0])
        colors = [colorFader("#cc2929", "#2929cc", i) for i in mixspace]
        # green 29cc2c
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        mpl.rcParams["font.size"] = 32
        fig, plot1 = plt.subplots(1, 1)
        fig.set_dpi(50)
        fig.set_figwidth(5)
        fig.set_figheight(4)

        self.tabs_frame = ctk.CTkTabview(self)
        self.tabs_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        #print(dataname)
        #print(len(np.shape(data[1][0])))
        if dataname == "Phasecomp":
            fig.gca().set_prop_cycle('color', colors)
            for i in range(len(data[1])):
                plot1.plot(data[0], data[1][i])
            plot1.legend(["Martensite", "Austenite", "Bainite", "Ferrite", "Pearlite"])
        elif dataname == "Matcomp":
            fig.gca().set_prop_cycle('color', colors)
            for i in range(len(data[1])):
                plot1.plot(data[0], data[1][i])
            plot1.legend(["C", "N", "Cr", "Ni", "Si"])
        elif len(np.shape(data[1])) == 3:
            #print("Testing")
            #print("Depth in xyz" + str(np.shape(data[1][:, 0, 0])))
            """ Data for ponts"""
            fig.gca().set_prop_cycle('color', colors)

            if len(data[1][0, :, 0]) == 6:
                direction = ["x", "y", "z", "yz", "xz", "xy"]
            elif len(data[1][0, :, 0]) == 3:
                direction = ["1", "2", "3"]
            else:
                direction = np.linspace(0, len(data[1][0, :, 0]), len(data[1][0, :, 0]))

            for i in range(len(data[1][0, :, 0])):
                # Do a stack to go one level deeper
                tab0 = self.tabs_frame.add(direction[i])
                tab0frame = resultTab(tab0, dataname+"_"+str(i), [data[0], data[1][:, i, :]], leg, xlbl, ylbl)
                tab0.rowconfigure(0, weight=1)
                tab0.columnconfigure(0, weight=1)
                tab0frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
            return

        elif len(np.shape(data[1][0])) == 2:
            fig.gca().set_prop_cycle('color', colors)
            plot1.plot(data[0], data[1])
            plot1.legend(leg)
        else:
            fig.gca().set_prop_cycle('color', colors)
            plot1.plot(data[0], np.transpose(data[1]))
            plot1.legend(leg)

        plot1.set_xlabel(xlbl)
        plot1.set_ylabel(ylbl)


        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")
class headerFrame(ctk.CTkFrame):
    def __init__(self, master, filename):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)
        self.filetextbox= ctk.CTkTextbox(self, width=400, height=25)
        self.filetextbox.insert("0.0", filename)
        self.filetextbox.grid(row=0, column=0, padx=20, pady=20, sticky="e")
        self.sidebar_button_1 = ctk.CTkButton(self, text="Continue")
        self.sidebar_button_1.grid(row=0, column=1, padx=20, pady=20, sticky="w")

class resultFrame(ctk.CTkFrame):
    """
    Frame for result information. The information is structured in a tab view with graphs for each dataset in xdmf file.

    The Tabs are automatically constructed from datasets. Composition

    If phase fractions are in the dataset another tab with the final timestep is constructed over the full geometry.

    """

    def __init__(self, master, filename):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.tabs_frame = ctk.CTkTabview(self)
        # self.tabs_frame.rowconfigure(0, weight=1)
        # self.tabs_frame.columnconfigure(0, weight=1)
        self.tabs_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")


        tabs = ["Austenite", "Bainite", "Martensite", "Stress", "vonMises", "Plastic strain", "Elastic strain", "Temperature", "Composition/C"]
        datanames = ["Austenite", "Bainite", "Martensite", "Stress", "vonMises", "Strain_pl", "Strain", "T", "Composition/C"]

        excl = ['JMAK_tau_Ferrite', 'JMAK_n_Ferrite', 'JMAK_tau_Bainite', 'JMAK_n_Bainite', 'JMAK_tau_Pearlite', 'JMAK_n_Pearlite', 'KM_Ms_Martensite', 'KM_b_Martensite']
        excl_comp = ["Composition/C", "Composition/N", 'Composition/Cr', 'Composition/Mn', 'Composition/Ni', 'Composition/Mo', 'Composition/Si']
        excl_phases = ["Austenite", "Bainite", "Ferrite", "Pearlite"]
        tabs = list(getnames_results(filename))
        tabs = [i for i in tabs if i not in excl and i not in excl_comp and i not in excl_phases]
        allpoints = read_results(filename, "nodes")

        radius = np.max(allpoints[:, 0])
        points = [[0, 0, 0],
             [radius / 2, 0, 0],
             [6 * radius / 10, 0, 0],
             [7 * radius / 10, 0, 0],
             [8 * radius / 10, 0, 0],
             [9 * radius / 10, 0, 0],
             [9.5 * radius / 10, 0, 0],
             [radius, 0, 0]]
        points_leg = np.array(points)[:, 0]*1000
        points_leg = points_leg.round(1)
        # Changing the point data to 2D if the dataset is 2D
        if np.shape(allpoints)[1] == 2:
            points = np.array(points)[:, 0:2]


        alldata_dict = read_results_all(filename, points)
        if "Martensite" in tabs:
            tab0 = self.tabs_frame.add("Phase composition")
            phases_data = list()
            phases = ["Martensite", "Austenite", "Bainite", "Ferrite", "Pearlite"]
            for phase in phases:
                phases_data.append(read_results_axis(filename, phase, -1))
            pcompdata = [read_results_axis(filename, "nodes")[:, 0], phases_data]
            tab0frame = resultTab(tab0, "Phasecomp", pcompdata, phases)
            tab0.rowconfigure(0, weight=1)
            tab0.columnconfigure(0, weight=1)
            tab0frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        tab0 = self.tabs_frame.add("Material composition")
        compositions = list()
        elements = ["C", "N", "Cr", "Ni", "Si"]
        for el in elements:
            compositions.append(read_results_axis(filename, "Composition/" + el))
        pcompdata = [read_results_axis(filename, "nodes")[:, 0], compositions]
        tab0frame = resultTab(tab0, "Matcomp", pcompdata, elements)
        tab0.rowconfigure(0, weight=1)
        tab0.columnconfigure(0, weight=1)
        tab0frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        for i in range(len(tabs)):
            if len(np.shape(alldata_dict[tabs[i]][1])) == 3:
                tmpdata = np.transpose(np.array(alldata_dict[tabs[i]][1]), (1, 2, 0))
            elif len(np.shape(alldata_dict[tabs[i]][1])) == 2:
                tmpdata = np.transpose(np.array(alldata_dict[tabs[i]][1]))
            tab0 = self.tabs_frame.add(tabs[i])
            tab0frame = resultTab(tab0, tabs[i], [alldata_dict[tabs[i]][0], tmpdata], points_leg)
            tab0.rowconfigure(0, weight=1)
            tab0.columnconfigure(0, weight=1)
            tab0frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        print("Displayed results")
        #self.tabs_frame.set("Composition/C")

class Result_MainApp(ctk.CTk):
    def __init__(self, filename):

        super().__init__()
        self.wm_iconbitmap("GUIfiles/Ballbearing.ico")
        self.geometry("1600x1000")
        self.title("Quenching of steel")

        # Setting weights
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)

        # Create header
        self.header_frame = headerFrame(self, filename)
        self.header_frame.grid(row=0, column=0, sticky="nsew")

        # Create sidebar
        self.sidebar_frame = resultFrame(self, filename)
        self.sidebar_frame.grid(row=1, column=0, sticky="nsew")


        # Adding button functionality
        self.header_frame.sidebar_button_1.configure(command=self.change_result)

    def change_result(self):
        """
        Changing the displayed results
        """
        filename = filedialog.askopenfilename(filetypes=(("xdmf files", "*.xdmf"),))
        self.header_frame.filetextbox.delete("0.0", "end")
        self.header_frame.filetextbox.insert("0.0", filename)
        newfilename = self.header_frame.filetextbox.get("0.0", ctk.END).strip("\n")
        print(newfilename)
        self.sidebar_frame.destroy()
        self.sidebar_frame = resultFrame(self, newfilename)
        self.sidebar_frame.grid(row=1, column=0, sticky="nsew")