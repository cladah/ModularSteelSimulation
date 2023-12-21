import tkinter as tk
from tkinter import ttk
class MainApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        #tk.Tk.iconbitmap(self, default="clienticon.ico")
        tk.Tk.wm_title(self, "Client")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        frame = StartPage(container, self)

        self.frames[StartPage] = frame

        frame.grid(row=0, column=0, sticky="nsew")


        self.show_frame(StartPage)
    def create_frame(self, cont):
        pass
        #self.frames[cont] = Page2(container, self)
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
class Page2(tk.Frame):
    pass
class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        from HelpFile import read_input
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        data = read_input()
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Start Page")
        label.pack(pady=10,padx=10)
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
        fig = Figure(figsize=(5, 4), dpi=100)
        plot1 = fig.add_subplot(111)
        plot1.plot(times, temps)
        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().pack()

def runguimodule(gui):
    import tkinter as tk
    #from tkinter import ttk
    from HelpFile import read_input
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                   NavigationToolbar2Tk)
    import threading
    import matplotlib.pyplot as plt
    data = read_input()

    #tabControl = tk.ttk.Notebook(gui)
    #tab1 = tk.ttk.Frame(tabControl)
    #tab2 = tk.ttk.Frame(tabControl)
    #tabControl.add(tab1, text='Input')
    #tabControl.pack(expand=1, fill="both")


    composition = tk.Label(gui, text="Composition of steel is:", pady=20)
    composition.pack()
    compstr = " ".join([i + "=" + str(data["Material"]["Composition"][i]) for i in data["Material"]["Composition"].keys()])
    composition = tk.Label(gui, text=compstr)
    composition.pack()

    text = tk.Label(gui, text="Temperature history:")
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
    canvas = FigureCanvasTkAgg(fig, master=gui)
    canvas.draw()
    canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    canvas.get_tk_widget().pack()


    return gui

def addgui(gui, type):
    from HelpFile import readdatastream, getaxisvalues, getTTTdata, read_input
    import numpy as np
    import tkinter as tk
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                   NavigationToolbar2Tk)


    if type == "Mesh":

        tabControl = tk.ttk.Notebook(gui)
        tab2 = tk.ttk.Frame(tabControl)
        tabControl.add(tab2, text='Input')
        pass
        # fig = Figure(figsize=(5, 4), dpi=100)
        # reader = vtk.vtkUnstructuredGridReader()
        # reader.SetFileName("Resultfiles/Mesh.vtk")
        # reader.Update()
        # data = reader.GetOutput()
        #
        # points = data.GetPoints()
        # npts = points.GetNumberOfPoints()
        # x = vtk_to_numpy(points.GetData())
        # triangles = vtk_to_numpy(data.GetCells().GetData())
        # ntri = triangles.size // 4  # number of cells
        # tri = np.take(triangles, [n for n in range(triangles.size) if n % 4 != 0]).reshape(ntri, 3)
        # plt.figure(figsize=(8, 8))
        # plt.triplot(x[:, 0], x[:, 1], tri)
        # plt.gca().set_aspect('equal')
        # plt.show()
    elif type == "Carbonitriding":
        wC = getaxisvalues("Composition/C")
        wN = getaxisvalues("Composition/N")
        xyz = getaxisvalues("nodes")
        fig = Figure(figsize=(5, 4), dpi=50)
        plot1 = fig.add_subplot(111)
        plot1.plot(np.array(xyz)[:,0], wC)
        plot1.plot(np.array(xyz)[:, 0], wN)
        canvas = FigureCanvasTkAgg(fig, master=gui)
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    elif type == "TTT":
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
        canvas = FigureCanvasTkAgg(fig, master=gui)
        canvas.draw()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        canvas.get_tk_widget().pack()


def _run():
    pass
