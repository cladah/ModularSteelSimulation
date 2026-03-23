import pathlib
import os
import json
import h5py
from sqlitedict import SqliteDict
import numpy as np
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt

def make_hashable_dict(d):
    return frozenset(sorted(d.items()))

def get_file_path(title="Select File", filetypes=(("XDMF files", "*.xdmf"), ("All files", "*.*"))):
    """Helper to handle file dialogs without cluttering the main logic."""
    root = tk.Tk()
    root.withdraw()  # Hide the main tkinter window
    path = filedialog.askopenfilename(title=title, filetypes=filetypes)
    root.destroy()
    return path

def read_old_input(filename):
    f = open(filename,'r')
    data = json.load(f)
    return data
def read_input():
    inputdata = read_geninput()
    for file in inputdata["Inputs"]:
        ftmp = open("Inputs/" + inputdata["InputDirectory"] + "/" + file + ".json", 'r')
        tmpdata = json.load(ftmp)
        inputdata.update(tmpdata)
        ftmp.close()
    return inputdata

def change_inputfile(filename):
    f = open(filename, 'r')
    data = json.load(f)
    f.close()

    f = open('Cachefiles/Input.json', 'w')
    json.dump(data, f, indent=2)
    f.close()

def reset_input(inputfile):
    f = open(inputfile, 'r')
    data = json.load(f)
    f.close()

    f = open('Cachefiles/Input.json', 'w')
    json.dump(data, f, indent=2)
    f.close()
def createinputcache():
    f = open('Cachefiles/InputCache.json', 'w')
    totinput = dict()
    ginput = read_geninput()
    totinput.update(ginput)
    for file in ginput["Inputs"]:
        minput = read_modinput("Inputs/"+ginput["InputDirectory"] + "/" + file + ".json")
        totinput.update(minput)
    json.dump(totinput, f, indent=2)
    f.close()

def createresultinput():
    ginput = read_geninput()
    filename = ginput["Datastream"]["Savedirect"].split('.')[0]
    f = open("Resultfiles/" + filename + ".json", "w")
    totinput = dict()
    totinput.update(ginput)
    for file in ginput["Inputs"]:
        minput = read_modinput("Inputs/"+ginput["InputDirectory"] + "/" + file + ".json")
        totinput.update(minput)
    json.dump(totinput, f, indent=2)
    f.close()

def adjustinputcache(model):
    f = open('Cachefiles/Input.json', 'r')
    indata = json.load(f)
    f.close()
    f = open('Cachefiles/InputCache.json', 'w')
    data = read_input()
    if model == 'Mesh':
        data['Geometry'] = indata['Geometry']
    elif model == 'Quenching':
        for x in ['Geometry', 'Material', 'Thermo', 'FEM', 'Programs']:
            data[x] = indata[x]
    else:
        raise KeyError('Adjust input cache not implemented')
    json.dump(data, f, indent=2)
    f.close()

def change_input(firstlvl, secondlvl, data):
    newdata = read_input()
    newdata[firstlvl][secondlvl] = data

    f = open('Cachefiles/Input.json', 'w')
    json.dump(newdata, f, indent=2)
    f.close()

def inputfile(filename):
    f = open(filename, 'r')
    newdata = json.load(f)
    f.close()
    f = open("Setup/ReferenceInput", 'r')
    reference = json.load(f)
    f.close()


    if newdata.keys() == reference.keys():
        f = open('Cachefiles/Input.json', 'w')
        json.dump(newdata, f, indent=2)
        f.close()
def checkruncondition(model):
    ginput = read_geninput()
    f = open('Cachefiles/InputCache.json', 'r')
    cachedata = json.load(f)
    f.close()

    modulelist = list(ginput["Modules"])
    # Check rerun criteria
    if ginput["RerunAll"] == True:
        return True

    for i in range(len(modulelist)):
        if ginput["Modules"][i] == model and ginput["Rerun"][i] == True:
            return True

    return False


def saveresult(filename, dataname,data):
    try:
        with h5py.File("Resultfiles/" + filename, "r+") as f:
            pass
    except:
        with h5py.File("Resultfiles/" + filename, "w") as f:
            pass

    with h5py.File("Resultfiles/"+filename, "r+") as f:
        try:
            del f[dataname]
        except:
            pass
        f.create_dataset(dataname, data=data)

def readresultfile(filename, dataname):
    try:
        with h5py.File("Resultfiles/"+filename, "r") as f:
            data = np.array(f.get(dataname))
        return data
    except:
        raise KeyError("Result \""+str(dataname)+"\" doesn't exist in result file")
def addTTTdb(data, comp_data, model_input, data_name):
    for key in comp_data:
        if key == "C" or key == "N":
            comp_data[key] = float(round(0.05 * round(comp_data[key]/0.05), 2))
        elif comp_data[key] != round(comp_data[key], 1):
            comp_data[key] = float(round(comp_data[key], 1))
    dbkey = str((sorted(comp_data.items()), sorted(model_input.items())))

    with SqliteDict("MaterialDatabase/database.db", autocommit=False) as db:
        entry = db.get(dbkey, {})

        if data_name in entry:
            print(f"Data for '{data_name}' already exists. Skipping.")
            return

        entry[data_name] = data
        db[dbkey] = entry

        db.commit()
        print(f"Added {data_name} to database for {comp_data}")

def getTTTdb(comp_data, model_input, data_name):
    for key in comp_data:
        if key == "C" or key == "N":
            comp_data[key] = float(round(0.05 * round(comp_data[key]/0.05), 2))
        elif comp_data[key] != round(comp_data[key], 1):
            comp_data[key] = (round(comp_data[key], 1))

    dbkey = str((sorted(comp_data.items()), sorted(model_input.items())))
    with SqliteDict("MaterialDatabase/database.db", flag='r') as db:
        # Check if the key exists AND if the specific type is inside it
        if dbkey in db:
            entry = db[dbkey]
            if data_name in entry:
                print(f"{data_name} exists for composition {comp_data}")
                return entry[data_name]

    print(f"{model_input}")
    raise KeyError(f"Database doesn't have information regaring {data_name} for composition {comp_data}")



    db = SqliteDict("MaterialDatabase/database.db", outer_stack=False)
    if dbkey in db:
        print("Data in database, checking type")
        if type in db[dbkey].keys():
            return db[dbkey][data_name]

    else:
        pass
    print(f"{model_input}")
    raise KeyError(f"Database doesn't have information regaring {data_name} for composition {comp_data}")

def checkTTTdb(comp_data, model_input, data_name):
    for key in comp_data:
        if key == "C" or key == "N":
            comp_data[key] = float(round(0.05 * round(comp_data[key]/0.05), 2))
        elif comp_data[key] != round(comp_data[key], 1):
            comp_data[key] = float(round(comp_data[key], 1))

    dbkey = str((sorted(comp_data.items()), sorted(model_input.items())))
    with SqliteDict("MaterialDatabase/database.db", flag='r') as db:
        # Check if the key exists AND if the specific type is inside it
        if dbkey in db:
            entry = db[dbkey]
            if data_name in entry:
                print(f"{data_name} exists for composition {comp_data}")
                return True
            else:
                return False
    return False
def checkDatabase():
    with SqliteDict("MaterialDatabase/database.db", flag='r') as db:
        # Check if the key exists AND if the specific type is inside it
        x = list()
        y = list()
        for dbkey in db:
            print(" ")

            safe_dict = {"np": np}

            data = eval(dbkey, {"__builtins__": {}}, safe_dict)

            first_float = data[0][0][1]

            x.append(first_float)
            y.append(db[dbkey]["Martensite_"]["start"][0][0])

        plt.plot(x,y,'o')
        plt.show()

def analyseTTTdatabase():

    TTTdata = SqliteDict("MaterialDatabase/database.db", tablename="TTTdata", outer_stack=False)
    print("Compositions in database")
    for key in TTTdata.keys():
        print(TTTdata[key]["TTTdata"].keys())
        #print(TTTdata[key]["Composition"])
        #print(TTTdata[key]["Modeldata"])
        #print(TTTdata[key]["TTTdata"])


def get_plotlbls(dataname):
    xlbls = ["Radius [mm]", "Time [s]"]
    ylbls = ["Weight [%]", "Stress [Pa]", "Strain [-]", "? [?]"]

    if "Composition" in dataname:
        xlbl = xlbls[0]
        ylbl = ylbls[0]
    elif dataname == "Stress":
        xlbl = xlbls[1]
        ylbl = ylbls[1]
    else:
        xlbl = xlbls[1]
        ylbl = ylbls[-1]

    return xlbl, ylbl

def read_geninput():
    """
    Reading general input files.

    Output:
        data - struct
    """
    fin = open("Inputs/input.json", 'r')
    inputdata = json.load(fin)
    fmain = open("Inputs/" + inputdata["InputDirectory"] + "/iMain.json", 'r')
    data = json.load(fmain)
    data.update(inputdata)
    fin.close()
    fmain.close()
    return data
def read_modinput(inputfile):
    f = open(inputfile, 'r')
    data = json.load(f)
    f.close()
    return data
def reset_output():
    ginput = read_geninput()
    from datetime import datetime
    with open("Resultfiles/output.txt", "w", encoding="utf-8") as file:
        file.write(datetime.now().strftime("%Y-%m-%d, %H:%M:%S") + "\nRunning Modular Steel Processing Simulation (MSPS)\n")
        file.write("Simulation software created by Clas Dahlin (2025)\n")

        file.write("---------------------------------------------------------------------\n\n")
        file.write("Modules used:\n")
        for module in ginput["Modules"]:
            file.write(str(module) + ": " + str(ginput["Programs"][module] + "\n"))
        file.write("---------------------------------------------------------------------\n\n")
        file.write("Data is saved using xdmf, with filename. " + str(ginput["Datastream"]["Savedirect"]) + "\n")
        file.write("---------------------------------------------------------------------\n\n")

def print_output(outputstr):
    # Opening and appending output string to textfile
    with open("Resultfiles/output.txt", "a", encoding="utf-8") as file:
        file.write(outputstr + "\n")