
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
        writer.write_points_cells(domain.points[:,:2], cells={"quad": domain.get_cells_type(list(domain.cells_dict.keys())[0])})
        writer.write_data(0, cell_data={})
    for element in ginput["Material"]["Composition"].keys():
        tmpelvalues = np.full(len(domain.points), ginput["Material"]["Composition"][element])
        adjustdatastream({"Composition_" + element: tmpelvalues}, "nodes")
    print("Created datastream file")
def adjustdatastream(data, datapos="nodes", t_data=0.0):
    dataname = list(data.keys())[0]

    if not isinstance(data, dict):
        raise KeyError("datatype for datastream storage should be a dict")
    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()
        n_steps = reader.num_steps
        point_data_list = []
        cell_data_list = []
        t_list = []
        print(f"File has {n_steps} timesteps:")
        for k in range(n_steps):
            t, point_data, cell_data = reader.read_data(k)
            point_data_list.append(point_data)
            cell_data_list.append(cell_data)
            t_list.append(t)
            print(f"  Step {k}: time={t}, point_data={list(point_data.keys())}")
    new = 0
    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        # write the mesh structure once
        writer.write_points_cells(points, cells)
        # Re-write old timesteps
        for k in range(n_steps):
            if t_data == t_list[k]:
                new = 1
                for dataname in data.keys():
                    point_data_list[k][dataname] = data[dataname]
            writer.write_data(t_list[k], point_data=point_data_list[k], cell_data=cell_data_list[k])
        if new == 0:
            writer.write_data(t_data, point_data=data)
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
        jsonfile = filename.split('.')[0] + ".json"
        pass
    elif filename.endswith(".h5"):
        jsonfile = filename.split('.')[0] + ".json"
        pass
    else:
        raise KeyError("Wrong format for saveing of datastream. Use .xdmf or .h5")
    try:
        with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
            points, cells = reader.read_points_cells()
            n_steps = reader.num_steps
            point_data_list = []
            cell_data_list = []
            t_list = []
            print(f"File has {n_steps} timesteps:")
            for k in range(n_steps):
                t, point_data, cell_data = reader.read_data(k)
                point_data_list.append(point_data)
                cell_data_list.append(cell_data)
                t_list.append(t)
                print(f"  Step {k}: time={t}, point_data={list(point_data.keys())}")
        org_dir = os.getcwd()
        os.chdir("Resultfiles")
        with meshio.xdmf.TimeSeriesWriter(filename) as writer:
            # write the mesh structure once
            writer.write_points_cells(points, cells)
            # Re-write old timesteps
            for k in range(n_steps):
                writer.write_data(t_list[k], point_data=point_data_list[k], cell_data=cell_data_list[k])
        os.chdir(org_dir)
        createresultinput()
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
    ginput = read_geninput()
    if ginput["Geometry"]["Type"]=="2D":
        node_y = readdatastream('nodes')[:, 1]
        indx = np.where(node_y == 0)
        data = readdatastream(dataname, time)[indx]
        x = readdatastream('nodes')[:, 0][indx]
        # Sorting data in order from x=min(x) to x=max(x)
        indx = np.argsort(x)
        data = np.array(data)[indx]
        if dataname == "nodes":
            data = data[:, 0]
    elif ginput["Geometry"]["Type"]=="4PointBend":
        node_x = readdatastream('nodes')[:, 0]
        indx = np.where(node_x == 0)
        data = readdatastream(dataname, time)[indx]
        node_y = readdatastream('nodes')[:, 1][indx]
        # Sorting data in order from x=min(x) to x=max(x)
        indx = np.argsort(node_y)
        data = np.array(data)[indx]
        if dataname == "nodes":
            data = data[:,1]
    else:
        raise KeyError("Geometry not implemented in getaxisvalues")
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

def readDatastreamMesh():
    with io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf", "r") as infile:
        # Read the mesh
        domain = infile.read_mesh()
    return domain
