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
    def __init__(self, master, dataname, data, points):
        super().__init__(master)

        filename = "Resultfiles/2024.xdmf"

        if "Composition" in dataname:
            xlbl = "Radius [mm]"
            ylbl = "Weight [%]"
        else:
            xlbl = "Time [s]"
            ylbl = "Stress [Pa]"



        xyz = read_results(filename, "nodes")
        radius = np.max(xyz[:, 0])
        names = getnames_results(filename)
        leg = []
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        mpl.rcParams["font.size"] = 32
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        for y in points:
            if len(np.shape(data[1])) != 3:
                plot1.plot(data[0], data[1])
            else:
                plot1.plot(data[0], np.array(data[1])[:,0])
                plot1.plot(data[0], np.array(data[1])[:,2])
                plot1.plot(data[0], np.array(data[1])[:,4])
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
        tabs = ["Stress", "vonMises", "Plastic strain", "Elastic strain", "Temperature", "Composition/C"]
        datanames = ["Stress", "vonMises", "Strain_pl", "Strain", "T", "Composition/C"]




        filename = "Resultfiles/2024.xdmf"
        xyz = read_results(filename, "nodes")
        radius = np.max(xyz[:, 0])
        points = [[radius / 2, 0],
             [6 * radius / 10, 0],
             [7 * radius / 10, 0],
             [8 * radius / 10, 0],
             [9 * radius / 10, 0],
             [radius, 0]]
        data_dict = read_results_all(filename, points)

        for i in range(len(tabs)):
            tab0 = self.tabs_frame.add(tabs[i])
            tab0frame = resultTab(tab0, datanames[i], data_dict[datanames[i]], points)
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