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

def reset_input():
    f = open('Cachefiles/Input_ref.json', 'r')
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
        if indata['Rerun'][m] == True:
            return True
        if m == model:
            break


    return False



    if indata["Programs"][model] != cachedata["Programs"][model]:
        return True

    if model == 'Meshing':
        if not os.path.isfile(pathlib.Path("Datastream.xdmf")):
            return True
        for x in ['Geometry']:
            if indata[x] != cachedata[x]:
                return True
    elif model == 'Carbonitriding':
        # pass
        for x in ['Geometry', 'Material', 'Thermo']:
            if indata[x] != cachedata[x]:
                return True
    elif model == 'TTT':
        for x in ['Geometry', 'Material', 'Thermo', 'Programs']:
            if indata[x] != cachedata[x]:
                return True
    elif model == 'Transformationmodels':
        for x in ['Geometry', 'Material', 'Thermo', 'Programs']:
            if indata[x] != cachedata[x]:
                return True
    elif model == 'Quenching':
        for x in ['Geometry', 'Material', 'Thermo', 'FEM', 'Programs']:
            if indata[x] != cachedata[x]:
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
def addTTTdata(compdata, data, type):
    if type not in ["TTTdata", "Modeldata"]:
        raise KeyError("Type can't be added to TTT database")
    from sqlitedict import SqliteDict
    TTTdata = SqliteDict("Resultfiles/database.db", tablename="TTTdata", outer_stack=False)
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
        if compdata[key] != round(compdata[key], 1):
            compdata[key] = round(compdata[key], 1)
    TTTdata = SqliteDict("Resultfiles/database.db", tablename="TTTdata", outer_stack=False)
    for key in TTTdata.keys():
        if compdata == TTTdata[key]["Composition"]:
            if type in TTTdata[key].keys():
                return TTTdata[key][type]
            else:
                raise KeyError(type + ' not in database for composition ' + str(compdata))
    print("Composition not in database")
def analyseTTTdatabase():

    TTTdata = SqliteDict("Resultfiles/database.db", tablename="TTTdata", outer_stack=False)
    print("Compositions in database")
    for key in TTTdata.keys():
        print(TTTdata[key]["Composition"])

def getTTTcompositions():
    roundingTTT = 1
    data = read_input()
    TTTcompositions = list()
    fullcomposition = dict()
    for element in data['Material']['Composition'].keys():
        fullcomposition[element] = getaxisvalues("Composition/" + element)
    composition = data['Material']['Composition']

    mesh = list()
    # Trying to create grid of compositions
    for element in data['Material']['Composition'].keys():
        if composition[element] != fullcomposition[element][-1]:
            if element == "C":
                tmplist = np.linspace(composition[element], fullcomposition[element][-1], 5)
            else:
                tmplist = np.linspace(composition[element], fullcomposition[element][-1], 2)
            tmplist = [round(elem, roundingTTT) for elem in tmplist]
            tmplist = list(set(tmplist)) # getting unique values
            tmplist.sort()
        else:
            tmplist = [composition[element]]
        mesh.append(tmplist)

    g = np.meshgrid(*mesh)
    positions = np.vstack(list(map(np.ravel, g)))
    for compnr in range(len(positions[0, :])):  # The number 0 here is correlated to the coal as it varies the most
        tmpcomp = dict()
        i = 0
        for element in data['Material']['Composition'].keys():
            tmpcomp[element] = round(positions[i, compnr], roundingTTT)
            i = i+1
        TTTcompositions.append(tmpcomp)
    return TTTcompositions

def setupSimulation():
    pass

