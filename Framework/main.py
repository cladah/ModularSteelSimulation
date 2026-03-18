import time
import tkinter
import tkinter as tk
import sys
import meshio
import json
from tkinter import filedialog
from pathlib import Path
from CGUImodule import MainApp
from Result_GUI import Result_MainApp, Compare_MainApp
from Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input, createinputcache, change_input, reset_input, analyseTTTdatabase, get_plotlbls, read_geninput, reset_output, get_file_path
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Diffusionmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule
from Modulefiles.MechanicalTest_file import MechTestModule
from Modulefiles.Testmod_file import TestModule
from Postprocessing.dataextraction import DatastreamPlotting, ResultPlotting, export_data, testread, interactive_datastream_plotter

def progressmonitor(tid, module):
    """
    :param tid:
    :param module:
    """
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)


def GUI():
    """
    Using Tkinter to run the GUI
    """
    print("Opening GUI window...")
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()

def Result_GUI_comp(filenames, dataname):
    """
            Using Tkinter to run the GUI
        """
    print("Opening result window...")
    ctk.set_appearance_mode("dark")
    app = Compare_MainApp(filenames, dataname)
    app.mainloop()

def Result_GUI_show(filename):
    """
        Using Tkinter to run the GUI
    """
    if filename == "":
        root = tkinter.Tk()
        root.withdraw()
        filename = tkinter.filedialog.askopenfilename(filetypes=(("xdmf files", "*.xdmf"),))
    print("Opening result window...")
    ctk.set_appearance_mode("dark")
    app = Result_MainApp(filename)
    app.mainloop()

def vtxfile_test():
    from adios2 import FileReader

    with FileReader("FeniCSx/beam_stress.bp") as s:
        # inspect variables
        vars = s.available_variables()
        for name, info in vars.items():
            print("variable_name: " + name, end=" ")
            for key, value in info.items():
                print("\t" + key + ": " + value, end=" ")
            print()
        print()
    pass


def confirm_simulation_config():
    # 1. Open Directory Selector
    initial_dir = Path.cwd() / "Inputs"
    if not initial_dir.exists():
        initial_dir = Path.cwd()

    root = tk.Tk()
    root.withdraw()  # Hide the annoying little white window
    selected_dir = filedialog.askdirectory(
        initialdir=str(initial_dir),
        title="Select Simulation Folder (within Inputs)"
    )
    root.destroy()

    if not selected_dir:
        print("Selection cancelled.")
        return None

    # 2. Locate and Load input.json
    config_path = Path(selected_dir) / "iMain.json"

    if not config_path.exists():
        print(f"Error: Could not find 'input.json' in {selected_dir}")
        return None
    try:
        with open(initial_dir / "input.json", "r") as f:
            config_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: {config_path.name} contains invalid JSON.")
        return None

    config_data["InputDirectory"] =  Path(selected_dir).name

    # 3. Display Info and Get Approval
    print("\n--- Configuration Found ---")
    print(json.dumps(config_data, indent=4))
    print("---------------------------\n")

    confirm = input(f"Does this look right for directory '{config_path.parent.name}'? (y/n): ").lower()

    if confirm == 'y':
        print("Config approved. Proceeding...")
        with open(initial_dir / "input.json", "w") as f:
            json.dump(config_data, f)
        return config_data
    else:
        print("Approval denied. Aborting.")
        return None


def setupSimulation():
    """
    Simulation setup from information in input.json

    :return:
    modules - Modules in order of exicution (list)
    """

    print("Setting up simulation...")
    reset_output()
    ginput = read_geninput()
    createdatastreamcache(ginput["Datastream"]["Cachedirect"])
    createinputcache()
    modules = list()
    for i in range(len(ginput["Modules"])):
        infile = ginput["InputDirectory"] + "/" + ginput["Inputs"][i]
        if ginput["Modules"][i] == "Meshing":
            modules.append(Meshingmodule(infile, i))
        elif ginput["Modules"][i] == "Test":
            modules.append(TestModule(infile, i))
        elif ginput["Modules"][i] == "Diffusion":
            modules.append(Diffusionmodule(infile, i))
        elif ginput["Modules"][i] == "TTTdiagram":
            modules.append(TTTdiagrammodule(infile, i))
        elif ginput["Modules"][i] == "TransformMod":
            modules.append(Transformationmodelmodule(infile, i))
        elif ginput["Modules"][i] == "Quenching":
            modules.append(Quenchingmodule(infile, i))
        elif ginput["Modules"][i] == "MechTest":
            modules.append(MechTestModule(infile, i))
        else:
            raise KeyError("Module input in iMain not supported")
    print("Simulation structure setup.")
    return modules

def main():
    menu = (
        "\n--- Simulation Manager ---\n"
        "1 - Run Simulation\n"
        "2 - Run using GUI\n"
        "3 - View Results\n"
        "4 - Export Data\n"
        "5 - Compare Files\n"
        "6 - Debug/Test Read\n"
        "0 - Exit\n"
        "Selection: "
    )
    while True:
        user_input = input(menu).strip()
        print(" ")
        if not user_input.isdigit():
            print("Invalid input. Please enter a number.")
            continue

        choice = int(user_input)
        match choice:
            case 1:
                print(">>> Starting normal input run...")
                config = confirm_simulation_config()
                if config:
                    ginput = read_geninput()
                    save_dir = ginput.get("Datastream", {}).get("Savedirect", "Default_Result")

                    modules = setupSimulation()
                    for module in modules:
                        module.run()
                        savedatastream(save_dir)

                    #Result_GUI_show(f"Resultfiles/{save_dir}")
                    break
            case 2:
                GUI()
            case 3:
                Result_GUI_show("")
            case 4:
                filename = get_file_path("Select XDMF to Export")
                if filename:
                    # Note: consider moving hardcoded strings to a config file
                    export_data(filename, ["Composition/C", "Martensite"], -1)

            case 5:
                file1 = get_file_path("Select First XDMF")
                file2 = get_file_path("Select Second XDMF")
                if file1 and file2:
                    ResultPlotting([file1, file2], "Composition/C")

            case 6:
                interactive_datastream_plotter()
                #testread("Datastream.xdmf")

            case 0:
                print("Exiting...")
                break

            case _:
                print("Option not recognized. Try again.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Goodbye!")
        sys.exit(0)
    print("End")
    sys.exit(0)
