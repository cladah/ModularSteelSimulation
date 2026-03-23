from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Diffusionmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule
from Modulefiles.MechanicalTest_file import MechTestModule
from Modulefiles.Optimization_file import OptimizationModule
from Modulefiles.Testmod_file import TestModule
from HelpFile import createinputcache, read_geninput, reset_output, get_file_path, checkDatabase
from Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
import json
from tkinter import filedialog
from pathlib import Path
import tkinter as tk
from scipy.optimize import minimize
from scipy.interpolate import interp1d

def update_config_json(newdirectory, filename="input.json"):
    file = Path.cwd() / "Inputs" / filename
    if file:
        with open(file, 'r+') as f:
            data = json.load(f)
            data['InputDirectory'] = newdirectory.name
            f.seek(0)
            json.dump(data, f, indent=4)
            f.truncate()
        return data
    return None

def run_simulation(config=False, inputdir=False):
    if inputdir:
        config = update_config_json(inputdir)
    if not config:
        config = confirm_simulation_config()
    if config:
        ginput = read_geninput()
        save_dir = ginput.get("Datastream", {}).get("Savedirect", "Default_Result")
        modules = setupSimulation()
        for module in modules:
            module.run()
            savedatastream(save_dir)
    return config


def optimization():
    confirm_simulation_config()
    ginput = read_geninput()
    infile = ginput["InputDirectory"] + "/" + ginput["Optimization"]
    OptMod = OptimizationModule(infile,0, run_simulation)
    OptMod.run()

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
def confirm_simulation_config_old():
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
        print(f"Error: Could not find 'iMain.json' in {selected_dir}")
        return None
    try:
        with open(initial_dir / "input.json", "r") as f:
            config_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: {config_path.name} contains invalid JSON.")
        return None

    config_data["InputDirectory"] = Path(selected_dir).name

    # Auto-update the directory name based on selection
    config_data["InputDirectory"] = Path(selected_dir).name

    # 3. Display Info and Get Approval Loop
    while True:
        print("\n----------------- Current Configuration ---------------------------")
        print(json.dumps(config_data, indent=4))
        print("---------------------------------------------------------------------\n")

        confirm = input("Does this look right? (y/n): ").lower()

        if confirm == 'y':
            print("Config approved. Saving and proceeding...")
            with open(initial_dir / "input.json", "w") as f:
                json.dump(config_data, f, indent=4)
            return config_data

        # 4. Editing Logic
        print("\n--- Edit Mode ---")
        print("Press Enter to keep the current value, or type a new value.")

        # Edit top-level keys
        for key in ["RerunAll", "InputDirectory"]:
            current_val = config_data[key]
            new_val = input(f"{key} [{current_val}]: ").strip()
            if new_val:
                # Basic type conversion for the boolean
                if key == "RerunAll":
                    config_data[key] = new_val.lower() in ['true', '1', 't', 'y']
                else:
                    config_data[key] = new_val

        # Edit nested Datastream keys
        for key in config_data["Datastream"]:
            current_val = config_data["Datastream"][key]
            new_val = input(f"Datastream -> {key} [{current_val}]: ").strip()
            if new_val:
                config_data["Datastream"][key] = new_val