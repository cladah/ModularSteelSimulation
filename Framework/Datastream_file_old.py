
import numpy as np
import meshio
from HelpFile import read_input, createinputcache, createresultinput, read_geninput
import os
import pathlib
import vtk
from dolfinx import io, fem
from mpi4py import MPI
import h5py
import ufl
import adios4dolfinx

class Datastream:
    def __init__(self, filename):
        self.filename = filename
        self.data = []

    def getstructure(self):
        return self.data

    def adddata(self, data, datapos="nodes", t_data=0):
        try:
            if type(data) == dict:
                pass
            else:
                raise TypeError
        except TypeError:
            raise KeyError("datatype for datastream storage should be a dict")

        with io.XDMFFile(MPI.COMM_WORLD,"Datastream.xdmf") as xdmf:
            xdmf.write_function(data,t_data)

    def getdata(self):
        pass

    def savetofile(self):
        pass

def createdatastream(domain):
    try:
        os.remove(os.getcwd() + "/Datastream.h5")
        os.remove(os.getcwd() + "/Datastream.xdmf")
    except FileNotFoundError:
        pass
    filename = os.getcwd() + "/Datastream.xdmf"

    ginput = read_geninput()
    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(domain.points[:,:2], cells={"quad": domain.get_cells_type("quad")})
        writer.write_data(0, cell_data={})
    # tmpelvalues = np.full(len(nodes), 0)
    for element in ginput["Material"]["Composition"].keys():
        tmpelvalues = np.full(len(domain.points), ginput["Material"]["Composition"][element])
        adjustdatastream({"Composition/" + element: tmpelvalues}, "nodes")

    print("Created datastream file")
def adjustdatastream(data, datapos="nodes", t_data=0.0):
    dataname = list(data.keys())[0]

    if not isinstance(data, dict):
        raise KeyError("datatype for datastream storage should be a dict")

    with io.XDMFFile(MPI.COMM_WORLD,"Datastream.xdmf", "r")as xdmf:
        domain = xdmf.read_mesh(name="Grid")
        #try:
        #xdmf.read_geometry_data(dataname)
        #except:
        #    t_data+=1.0

    if len(np.shape(data[dataname])) > 1:
        V = fem.functionspace(domain,("CG",2, (np.shape(data[dataname])[1],)))
    else:
        V = fem.functionspace(domain, ("CG", 2))
    u = fem.Function(V)
    u.name = list(data.keys())[0]
    u.x.array[:] = data[u.name].flatten()

    xdmf = io.XDMFFile(MPI.COMM_WORLD, 'Datastream.xdmf', 'a')
    xdmf.write_function(u, t_data)
    xdmf.close()

    return
def getnamesdatastream():
    with io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf") as xdmf:
        mesh = xdmf.read_mesh()
        with h5py.File("Datastream.h5", "r") as f:
            for k in f.keys():
                if k == "Function":
                    for k in f["Function"]:
                        print("  ", k)

    return

def readdatastream(dataname, time=0, all_t=0):
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
                if t == time and all_t == 0:
                    return point_data[dataname]
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)
            if all_t == 1:
                print("Getting all time values")
                return t_list, pd_list
            if time == -1:
                maxindx = t_list.index(np.max(t_list))
                return pd_list[maxindx][dataname]
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))
    return


def savedatastream(filename):
    if filename is None or filename == "":
        print("Datastream not saved. File in main folder.")
        return

    if filename.endswith(".xdmf"):
        pass
    elif filename.endswith(".h5"):
        pass
    else:
        raise KeyError("Wrong format for saveing of datastream. Use .xdmf or .h5")

    try:
        with io.XDMFFile(MPI.COMM_WORLD,"Datastream.xdmf", "r") as xdmf:
            domain = xdmf.read_mesh(name="Grid")
        V = fem.functionspace(domain, ("CG", 2))
        u = fem.Function(V)
        createinputcache()

        print("Saved datastream to " + filename)
    except:
        raise KeyError("No datastream to save, something went wrong")
    #try:
    #    os.remove("Datastream.h5")
    #    os.remove("Datastream.xdmf")
    #except:
    #    print("Couldn't remove old datastream")

def createdatastreamcache(filename=None):
    try:
        os.remove("Datastream_Cache.h5")
        os.remove("Datastream_Cache.xdmf")
        print("Old cache datastream removed")
    except FileNotFoundError:
        pass

    try:
        if filename is None or filename == "":
            with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
                points, cells = reader.read_points_cells()
                pd_list, cd_list, t_list = list(), list(), list()
                for k in range(reader.num_steps):
                    t, point_data, cell_data = reader.read_data(k)
                    t_list.append(t)
                    pd_list.append(point_data)
                    cd_list.append(cell_data)
            print("Using Datastream.xdmf as cache file")
        else:
            if filename.split(".")[1] in ["h5", "xdmf"]:
                pass
            else:
                raise FileNotFoundError
            filename = filename.split(".")[0]
            print("Using " + filename + ".xdmf as cache file")
            with meshio.xdmf.TimeSeriesReader(pathlib.Path(os.getcwd() + "/" + filename + ".xdmf")) as reader:
                points, cells = reader.read_points_cells()
                pd_list, cd_list, t_list = list(), list(), list()
                for k in range(reader.num_steps):
                    t, point_data, cell_data = reader.read_data(k)
                    t_list.append(t)
                    pd_list.append(point_data)
                    cd_list.append(cell_data)
        with meshio.xdmf.TimeSeriesWriter("Datastream_Cache.xdmf") as writer:
            writer.write_points_cells(points, cells)
            for i in range(len(t_list)):
                writer.write_data(t=t_list[i], point_data=pd_list[i], cell_data=cd_list[i])
        print("Datastream cache created\n")

    except meshio._exceptions.ReadError:
        print("No datastream caches, running empty simulation")
    except FileNotFoundError:
        print("Cache FileNotFound")
        print("No datastream caches, running empty simulation")


def removedatastreamcache():
    try:
        os.remove(os.getcwd() + "/Cachefiles/Datastream_Cache.h5")
        os.remove(os.getcwd() + "/Cachefiles/Datastream_Cache.xdmf")
    except:
        print("No datastream caches")


def readdatastreamcache(dataname, time=0):
    try:
        with meshio.xdmf.TimeSeriesReader("Datastream_Cache.xdmf") as reader:
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
    except:
        raise KeyError("Datastream " + str(dataname) + " doesn't exist in datastream file. Data that exist is "
                   + str(pd_list[0].keys()))

def getaxisvalues(dataname, time=0):
    node_y = readdatastream('nodes')[:, 1]
    indx = np.where(node_y == 0)
    data = readdatastream(dataname, time)[indx]
    x = readdatastream('nodes')[:, 0][indx]

    # Sorting data in order from x=min(x) to x=max(x)
    indx = np.argsort(x)
    data = np.array(data)[indx]
    return data
def gethistoryvalues(dataname, position):
    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()
        if dataname == "nodes":
            return points
        elif dataname == "elements":
            return cells
        pd_list, t_list = list(), list()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            if dataname in point_data.keys():
                t_list.append(t)
                pd_list.append(point_data[dataname][position])
    return t_list, pd_list
def resetdatastream():
    try:
        os.remove("Datastream.h5")
        os.remove("Datastream.xdmf")
        print("Old datastream removed\n")
    except FileNotFoundError:
        pass