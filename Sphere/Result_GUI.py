import customtkinter as ctk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
import numpy as np
import matplotlib as mpl

from ResultReading import read_results_history, read_results, getnames_results, read_results_all

class resultTab(ctk.CTkFrame):
    def __init__(self, master, dataname, data, points, allpoints):
        super().__init__(master)

        filename = "Resultfiles/2024.xdmf"

        datanames = ["Austenite", "Bainite", "Martensite", "Stress", "vonMises", "Strain_pl", "Strain", "T"]
        ynames = ["Phase fraction [-]", "Phase fraction [-]", "Phase fraction [-]", "Stress [Pa]", "von-Mises stress [Pa]",
             "Plastic strain [-]", "Elastic strain [-]", "Temperature [K]"]

        if "Composition" in dataname:
            ylbl = "Weight fraction [-]"
            xlbl = "Radius [mm]"
            leg = []
        elif dataname in datanames:
            ylbl = ynames[datanames.index(dataname)]
            xlbl = "Time [s]"
            leg = [str(round(i, 3)) for i in np.array(points)[:, 0]]
        else:
            ylbl = "?"
            xlbl = "Time [s]"
            leg = [str(round(i, 3)) for i in np.array(points)[:, 0]]


        pltcolors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        pltcolors = pltcolors[0:len(points)]


        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        mpl.rcParams["font.size"] = 32
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        if len(np.shape(data[1])) != 3:
            plot1.plot(data[0], data[1])
        else:
            plt.gca().set_prop_cycle(None)
            plot1.plot(data[0], np.array(data[1])[:, 0])
            plt.gca().set_prop_cycle(None)
            #plot1.plot(data[0], np.array(data[1])[:, 2], "-")
            plt.gca().set_prop_cycle(None)
            #plot1.plot(data[0], np.array(data[1])[:, 4], "--")
        plot1.set_xlabel(xlbl)
        plot1.set_ylabel(ylbl)
        plot1.legend(leg)

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(canvas, self, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky="nsew")
        toolbar.update()
        canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")
class headerFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.sidebar_button_1 = ctk.CTkButton(self, text="Continue")
        self.sidebar_button_1.grid(row=0, column=0, padx=20, pady=20)

class resultFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.tabs_frame = ctk.CTkTabview(self)
        # self.tabs_frame.rowconfigure(0, weight=1)
        # self.tabs_frame.columnconfigure(0, weight=1)
        self.tabs_frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")
        tabs = ["Austenite", "Bainite", "Martensite", "Stress", "vonMises", "Plastic strain", "Elastic strain", "Temperature", "Composition/C"]
        datanames = ["Austenite", "Bainite", "Martensite", "Stress", "vonMises", "Strain_pl", "Strain", "T", "Composition/C"]
        y = ["Phase fraction [-]", "Phase fraction [-]", "Phase fraction [-]", "Stress [Pa]", "vonMises stress [Pa]", "Plastic strain [-]", "Elastic strain [-]", "Temperature [K]",
                     "Weight fraction [-]"]



        filename = "Resultfiles/2024.xdmf"
        allpoints = read_results(filename, "nodes")
        radius = np.max(allpoints[:, 0])
        points = [[radius / 2, 0],
             [6 * radius / 10, 0],
             [7 * radius / 10, 0],
             [8 * radius / 10, 0],
             [9 * radius / 10, 0],
             [radius, 0]]
        data_dict = read_results_all(filename, points)
        for i in range(len(tabs)):

            if "Composition" in datanames[i]:
                tab0 = self.tabs_frame.add("Composition")
                tab0frame = resultTab(tab0, "Composition", data_dict["Composition/C"], points, allpoints)
            else:
                tab0 = self.tabs_frame.add(tabs[i])
                tab0frame = resultTab(tab0, datanames[i], data_dict[datanames[i]], points, allpoints)
            tab0.rowconfigure(0, weight=1)
            tab0.columnconfigure(0, weight=1)
            tab0frame.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        #self.tabs_frame.set("Composition/C")



class Result_MainApp(ctk.CTk):
    def __init__(self):

        super().__init__()
        self.wm_iconbitmap("GUIfiles/Ballbearing.ico")
        self.geometry("1600x1000")
        self.title("Quenching of steel")

        # Setting weights
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(0, weight=0)

        # Create header
        self.header_frame = headerFrame(self)
        self.header_frame.grid(row=0, column=0, sticky="nsew")

        # Create sidebar
        self.sidebar_frame = resultFrame(self)
        self.sidebar_frame.grid(row=1, column=0, sticky="nsew")


        # Adding button functionality
        self.header_frame.sidebar_button_1.configure(command=self.next_result)

    def next_result(self):
        pass