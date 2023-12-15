
def runguimodule():
    import tkinter as tk
    from HelpFile import read_input
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                   NavigationToolbar2Tk)
    import threading
    import matplotlib.pyplot as plt
    data = read_input()
    gui = tk.Tk()
    gui.geometry("500x600")
    gui.title("Quenching of steel")
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
    tempertemp = 400 #data["Thermo"]["tempertemp"]
    roomtemp = data["Thermo"]["quenchtemp"]-273.15
    holdCN = data["Thermo"]["CNtime"]*60
    holdquench = holdCN + 30
    holdtemper = holdquench + 60
    holdend = holdtemper + 60
    times = [0, 0, holdCN,holdCN+10,holdquench,holdquench + 10, holdtemper,holdtemper + 10, holdend]
    temps = [roomtemp, starttemp, starttemp,quenchtemp,quenchtemp,tempertemp,tempertemp,roomtemp,roomtemp]
    fig = Figure(figsize=(5, 4), dpi=100)
    plot1 = fig.add_subplot(111)
    plot1.plot(times,temps)
    canvas = FigureCanvasTkAgg(fig, master=gui)
    canvas.draw()
    canvas.get_tk_widget().pack()



    return gui

def addgui(gui, type):
    import tkinter as tk
    button = tk.Button(master=gui, text="Quit", command=_run)
    button.pack(side=tk.RIGHT)
    gui.update()
    #if type=="carbonitriding":


def _run():
    pass