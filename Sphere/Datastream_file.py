import numpy as np
import meshio
from HelpFile import read_input, createinputcache
import os
import pathlib

def createdatastream():
    try:
        os.remove(os.getcwd() + "\Resultfiles\Datastream.h5")
        os.remove(os.getcwd() + "\Resultfiles\Datastream.xdmf")
    except FileNotFoundError:
        pass


def adjustdatastream(dataname, data, datapos="nodes", time=0):
    # Adding data to xdmf file
    if datapos == "nodes":
        meshstream = meshio.read("Resultfiles/Datastream.xdmf")
        meshstream.point_data[dataname] = data
        meshio.write("Resultfiles/Datastream.xdmf", meshstream)
    elif datapos == "elements":
        meshstream = meshio.read("Resultfiles/Datastream.xdmf")
        meshstream.cell_data[dataname] = data
        meshio.write("Resultfiles/Datastream.xdmf", meshstream)
    else:
        raise KeyError("datastream missing nodes or elements")


def readdatastream(dataname, time=0):
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
    if filename is None or filename == "":
        print("Datastream not saved. File in result folder.")
        return
    try:
        meshstream = meshio.read(os.getcwd() + "\Resultfiles\Datastream.xdmf")
        meshio.write(filename, meshstream)
        createinputcache()
        print("Saved datastream to " + filename)
    except meshio._exceptions.ReadError:
        raise KeyError("No datastream to save, something went wrong")


def createdatastreamcache(filename=None):
    try:
        os.remove(pathlib.Path("Cachefiles\\Datastream.h5"))
        os.remove(pathlib.Path("Cachefiles\\Datastream.xdmf"))
    except FileNotFoundError:
        pass

    try:
        if filename is None or filename == "":
            meshdata = meshio.read(pathlib.Path(os.getcwd() + "\\Resultfiles\\Datastream.xdmf"))
            meshio.write("Cachefiles\\Datastream.xdmf", meshdata)
        else:
            if filename.split(".")[1] in ["h5", "xdmf"]:
                pass
            else:
                raise FileNotFoundError
            filename = filename.split(".")[0]
            meshdata = meshio.read(pathlib.Path(os.getcwd() + "\\" + filename + ".xdmf"))
            meshio.write("Cachefiles\\Datastream.xdmf", meshdata)
        print("Datastream cache created\n")

    except meshio._exceptions.ReadError:
        print("No datastream caches, running empty simulation")
    except FileNotFoundError:
        print("Cache FileNotFound")
        print("No datastream caches, running empty simulation")


def removedatastreamcache():
    try:
        os.remove(os.getcwd() + "\Cachefiles\Datastream.h5")
        os.remove(os.getcwd() + "\Cachefiles\Datastream.xdmf")
    except:
        print("No datastream caches")


def readdatastreamcache(dataname, t=0):
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