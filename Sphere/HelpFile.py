import numpy as np
import meshio
import sqlite3
def read_input():
    import json
    f = open('Cachefiles/Input.json', 'r')
    data = json.load(f)
    f.close()
    return data
def createinputcache():
    import json
    f = open('Cachefiles/InputCache.json', 'w')
    data = read_input()
    json.dump(data, f, indent=2)
    f.close()

def adjustinputcache(model):
    import json
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

def checkinput(model):
    import json
    import numpy as np
    f = open('Cachefiles/Input.json', 'r')
    indata = json.load(f)
    f.close()
    f = open('Cachefiles/InputCache.json', 'r')
    cachedata = json.load(f)
    f.close()

    modellist = list(indata["Rerun"].keys())
    # create a stacked rerun criteria?

    # Check rerun criteria
    if indata['Rerun'][model] == 1:
        return False
    if model == 'Mesh':
        for x in ['Geometry']:
            if indata[x] != cachedata[x]:
                return False
    elif model == 'Carbonitriding':
        for x in ['Geometry', 'Material', 'Thermo', 'Programs']:
            if indata[x] != cachedata[x]:
                return False
    elif model == 'TTT':
        for x in ['Geometry', 'Material', 'Thermo', 'Programs']:
            if indata[x] != cachedata[x]:
                return False
    elif model == 'TTTfit':
        for x in ['Geometry', 'Material', 'Thermo', 'Programs']:
            if indata[x] != cachedata[x]:
                return False
    elif model == 'Quenching':
        for x in ['Geometry', 'Material', 'Thermo', 'FEM', 'Programs']:
            if indata[x] != cachedata[x]:
                return False
    return True

def createdatastream():
    data = read_input()
    # Geometry
    geometry = np.arange(data['Geometry']['nodes'])/(data['Geometry']['nodes']-1)
    # Composition / x-points
    comp0 = data['Material']['Composition']
    comp0list = list()
    for key in comp0.keys():
        comp0list.append(comp0[key])
    composition = [comp0list for i in range(data['Geometry']['nodes'])]
    # Temperature
    temperature = data['Geometry']['nodes']*[data['Thermo']['CNtemp']]
    # Phase fractions
    phases = [[1, 0, 0, 0] for i in range(data['Geometry']['nodes'])]
    # Empty parameter lists
    displacement = data['Geometry']['nodes'] * [0]
    TTT = data['Geometry']['nodes'] * [0]
    eps_pl = data['Geometry']['nodes'] * [0]
    eps_tr = data['Geometry']['nodes'] * [0]

    import h5py
    cachename = "Resultfiles/Datastream.hdf5"
    with h5py.File(cachename, "w") as f:
        f.create_dataset("geometry", data=geometry)
        f.create_dataset("composition", data=composition)
        f.create_dataset("temperature", data=temperature)
        f.create_dataset("phases", data=phases)
        f.create_dataset("TTT", data=TTT)
        f.create_dataset("displacement", data=displacement)
        f.create_dataset("plasticstrain", data=eps_pl)
        f.create_dataset("tripstrain", data=eps_tr)

def adjustdatastream(dataname,data,type):
    # Adding data to xdmf file
    if type == "nodes":
        meshstream = meshio.read("Resultfiles/Datastream.xdmf")
        meshstream.point_data[dataname] = data
        meshio.write("Resultfiles/Datastream.xdmf", meshstream)
    elif type == "elements":
        meshstream = meshio.read("Resultfiles/Datastream.xdmf")
        meshstream.cell_data[dataname] = data
        meshio.write("Resultfiles/Datastream.xdmf", meshstream)
    else:
        raise KeyError("datastream missing nodes or elements")
def readdatastream(dataname):
    import h5py
    try:
        meshstream = meshio.read("Resultfiles/Datastream.xdmf")
        if dataname in meshstream.point_data.keys():
            data = meshstream.point_data[dataname]
        elif dataname in meshstream.cell_data.keys():
            data = meshstream.cell_data[dataname]
        elif dataname == "nodes":
            data = meshstream.points
        elif dataname == "elements":
            data = meshstream.cells
        else:
            raise KeyError()
        return data
    except:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file")

def createdatastreamcache():
    import shutil
    import os
    try:
        shutil.copy("Resultfiles/Datastream.h5", "Cachefiles/Datastream.h5")
        shutil.copy("Resultfiles/Datastream.xdmf", "Cachefiles/Datastream.xdmf")
        os.remove("Resultfiles/Datastream.h5")
        os.remove("Resultfiles/Datastream.xdmf")
    except:
        print("No datastream caches")
def removedatastreamcache():
    import os
    try:
        os.remove("Cachefiles/Datastream.h5")
        os.remove("Cachefiles/Datastream.xdmf")
    except:
        print("No datastream caches")

def readdatastreamcache(dataname):
    import h5py
    try:
        meshstream = meshio.read("Cachefiles/Datastream.xdmf")
        if dataname in meshstream.point_data.keys():
            data = meshstream.point_data[dataname]
        elif dataname in meshstream.cell_data.keys():
            data = meshstream.cell_data[dataname]
        elif dataname == "nodes":
            data = meshstream.points
        elif dataname == "elements":
            data = meshstream.cells
        else:
            raise KeyError()
        return data
    except:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in cache file")
def getaxisvalues(dataname):
    node_y = readdatastream('nodes')[:, 1]
    indx = np.where(node_y == 0)
    y = readdatastream(dataname)[indx]
    x = readdatastream('nodes')[:, 0][indx]
    indx = np.argsort(x)
    y = np.array(y)[indx]
    return y

def saveresult(filename, dataname,data):
    import h5py
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
    import h5py
    import numpy as np
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
            if type in TTTdata[oldkey].keys():
                print(type + " for that composition is already stored in database")
                return
            print(type + " added to database for " + str(compdata))
            tmpdict = TTTdata[oldkey]
            tmpdict[type] = data
            TTTdata[oldkey] = tmpdict
            TTTdata.commit()
            TTTdata.close()
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
    from sqlitedict import SqliteDict
    TTTdata = SqliteDict("Resultfiles/database.db", tablename="TTTdata", outer_stack=False)
    for key in TTTdata.keys():
        if compdata == TTTdata[key]["Composition"]:
            if type in TTTdata[key].keys():
                return TTTdata[key][type]
            else:
                raise KeyError(type + ' not in database for composition ' + str(compdata))
    print("Composition not in database")
def analyseTTTdatabase():
    from sqlitedict import SqliteDict
    TTTdata = SqliteDict("Resultfiles/database.db", tablename="TTTdata", outer_stack=False)
    print("Compositions in database")
    for key in TTTdata.keys():
        print(TTTdata[key]["Composition"])