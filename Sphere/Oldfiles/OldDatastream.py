
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
def readdatastream(dataname, time=0):
    import h5py
    meshstream = meshio.read("Resultfiles/Datastream.xdmf")
    try:
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
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is " + str(*list(meshstream.point_data.keys())) + " and " + str(*list(meshstream.cell_data.keys())))

def savedatastream(filename):
    meshstream = meshio.read("Resultfiles/Datastream.xdmf")
    meshio.write(filename, meshstream)
    createinputcache()

def createdatastreamcache(filename=None):
    import shutil
    import os
    try:
        print(filename)
        if filename == None:
            meshdata = meshio.read("Resultfiles/Datastream.xdmf")
            meshio.write("Cachefiles/Datastream.xdmf", meshdata)
            os.remove("Resultfiles/Datastream.h5")
            os.remove("Resultfiles/Datastream.xdmf")
        else:
            # Add check if filename have xdmf or h5 extension
            filename = filename.split(".")[0]
            meshdata = meshio.read(filename + ".xdmf")
            meshio.write("Cachefiles/Datastream.xdmf", meshdata)
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
    meshstream = meshio.read("Cachefiles/Datastream.xdmf")
    try:
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
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in cache file. Data that exist is " + str(*list(meshstream.point_data.keys())) + " and " + str(*list(meshstream.cell_data.keys())))
def getaxisvalues(dataname, time=0):
    node_y = readdatastream('nodes')[:, 1]
    indx = np.where(node_y == 0)
    y = readdatastream(dataname)[indx]
    x = readdatastream('nodes')[:, 0][indx]
    indx = np.argsort(x)
    y = np.array(y)[indx]
    return y