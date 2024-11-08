import pathlib
import os
import json
import h5py
from sqlitedict import SqliteDict
import numpy as np


def read_input():
    f = open('Cachefiles/Input.json', 'r')
    data = json.load(f)
    f.close()
    return data

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
    data = read_input()
    json.dump(data, f, indent=2)
    f.close()

def createresultinput(filename):
    f = open("Resultfiles/" + filename + ".json", "w")
    data = read_input()
    json.dump(data, f, indent=2)
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
    f = open('Cachefiles/Input.json', 'r')
    indata = json.load(f)
    f.close()
    f = open('Cachefiles/InputCache.json', 'r')
    cachedata = json.load(f)
    f.close()

    modellist = list(indata["Rerun"].keys())
    # Check rerun criteria
    if indata["Rerun"]["All"] == True:
        return True

    for m in modellist:
        if indata['Rerun'][model] == True:
            return True
        #if m == model:
        #    break

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
def addTTTdata(compdata, data, type):
    if type not in ["TTTdata", "Modeldata"]:
        raise KeyError("Type can't be added to TTT database")
    from sqlitedict import SqliteDict
    TTTdata = SqliteDict("MaterialDatabase/database.db", tablename="TTTdata", outer_stack=False)
    for oldkey in TTTdata.keys():
        if compdata == TTTdata[oldkey]["Composition"]:
            # if type in TTTdata[oldkey].keys():
            #     print(type + " for that composition is already stored in database")
            #     return

            tmpdict = TTTdata[oldkey]
            tmpdict[type] = data
            TTTdata[oldkey] = tmpdict
            TTTdata.commit()
            TTTdata.close()
            print(type + " added to database for " + str(compdata))
            return
    print("Composition not in TTT database")
    TTTdata[str(len(TTTdata) + 1)] = {type:data,"Composition":compdata}
    TTTdata.commit()
    TTTdata.close()
    print(type + " is now stored in TTT database for " + str(compdata))
def getTTTdata(compdata, type):
    # Making sure the composition has the correct round value
    for key in compdata:
        if key == "C" or key == "N":
            compdata[key] = round(0.05 * round(compdata[key]/0.05), 2)
        elif compdata[key] != round(compdata[key], 1):
            compdata[key] = round(compdata[key], 1)
    TTTdata = SqliteDict("MaterialDatabase/database.db", tablename="TTTdata", outer_stack=False)
    for key in TTTdata.keys():
        if compdata == TTTdata[key]["Composition"]:
            if type in TTTdata[key].keys():
                return TTTdata[key][type]
            else:
                raise KeyError(type + ' not in database for composition ' + str(compdata))
    print("Composition not in database")
def analyseTTTdatabase():

    TTTdata = SqliteDict("MaterialDatabase/database.db", tablename="TTTdata", outer_stack=False)
    print("Compositions in database")
    for key in TTTdata.keys():
        print(TTTdata[key]["TTTdata"].keys())
        #print(TTTdata[key]["Composition"])
        #print(TTTdata[key]["Modeldata"])
        #print(TTTdata[key]["TTTdata"])

def setupSimulation():
    pass


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
