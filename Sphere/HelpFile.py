import numpy as np
import meshio
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
    # Check rerun criteria
    if indata['Rerun'][model] == 1:
        return False
    if model == 'Mesh':
        for x in ['Geometry']:
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



def saveresult(filename, dataname,data):
    import h5py
    with h5py.File("Resultfiles/"+filename, "r+") as f:
        try:
            del f[dataname]
        except:
            pass
        f.create_dataset(dataname, data=data)

def readresultfile(filename, dataname):
    import h5py
    try:
        with h5py.File("Resultfiles/"+filename, "r") as f:
            data = np.array(f.get(dataname))
        return data
    except:
        raise KeyError("Result "+str(dataname)+" doesn't exist in result file")