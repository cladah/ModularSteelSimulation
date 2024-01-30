import numpy as np
import meshio
from HelpFile import read_input, createinputcache
import os
import pathlib
#import xdmf

def createdatastream():
    try:
        os.remove(os.getcwd() + "/Datastream.h5")
        os.remove(os.getcwd() + "/Datastream.xdmf")
    except FileNotFoundError:
        pass


def adjustdatastream(dataname, data, datapos="nodes", t_data=0):
    # Adding data to xdmf file

    # meshstream = meshio.read("Datastream.xdmf")
    # if time is None:
    #     print("2")
    #     if datapos == "nodes":
    #         meshstream.point_data[dataname] = data
    #         meshio.write("Datastream.xdmf", meshstream)
    #     elif datapos == "elements":
    #         meshstream.cell_data[dataname] = data
    #         meshio.write("Datastream.xdmf", meshstream)
    #     else:
    #         raise KeyError("datastream missing nodes or elements")
    # else:
    try:
        with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
            points, cells = reader.read_points_cells()
            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
    except meshio._exceptions.ReadError:
        raise KeyError("No meshgrid for timeseries")
    if t_data not in t_list:
        t_list.append(t_data)
        pd_list.append({dataname: data})
        cd_list.append({})
    else:
        indx = t_list.index(t_data)
        pd_list[indx][dataname] = data
    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(points, cells)
        for i in range(len(t_list)):
            writer.write_data(t=t_list[i], point_data=pd_list[i], cell_data=cd_list[i])
    #print("Added " + str(dataname) + " to datastream")
    return

def readdatastream(dataname, time=0):
    try:
        with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
            points, cells = reader.read_points_cells()
            if dataname == "nodes":
                return points
            elif dataname == "elements":
                return cells
            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                if t == time:
                    return point_data[dataname]
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))
    return



    meshstream = meshio.read("Datastream.xdmf")
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
        meshstream = meshio.read(os.getcwd() + "/Datastream.xdmf")
        meshio.write(filename, meshstream)
        createinputcache()
        print("Saved datastream to " + filename)
    except meshio._exceptions.ReadError:
        raise KeyError("No datastream to save, something went wrong")


def createdatastreamcache(filename=None):
    try:
        os.remove("/Cachefiles/Datastream.h5")
        os.remove("/Cachefiles/Datastream.xdmf")
    except FileNotFoundError:
        pass

    try:
        if filename is None or filename == "":
            meshdata = meshio.read(pathlib.Path(os.getcwd() + "/Datastream.xdmf"))
            meshio.write("Cachefiles/Datastream.xdmf", meshdata)
        else:
            if filename.split(".")[1] in ["h5", "xdmf"]:
                pass
            else:
                raise FileNotFoundError
            filename = filename.split(".")[0]
            meshdata = meshio.read(pathlib.Path(os.getcwd() + "/" + filename + ".xdmf"))
            meshio.write("Cachefiles/Datastream.xdmf", meshdata)
        print("Datastream cache created\n")

    except meshio._exceptions.ReadError:
        print("No datastream caches, running empty simulation")
    except FileNotFoundError:
        print("Cache FileNotFound")
        print("No datastream caches, running empty simulation")


def removedatastreamcache():
    try:
        os.remove(os.getcwd() + "/Cachefiles/Datastream.h5")
        os.remove(os.getcwd() + "/Cachefiles/Datastream.xdmf")
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
    data = readdatastream(dataname)[indx]
    x = readdatastream('nodes')[:, 0][indx]
    indx = np.argsort(x)
    data = np.array(data)[indx]
    return data