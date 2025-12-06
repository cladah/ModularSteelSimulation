
import numpy as np
import meshio
from HelpFile import read_input, createinputcache, createresultinput, read_geninput
import os
import pathlib
import vtk
from dolfinx import io, fem, mesh as dmesh
from mpi4py import MPI
import h5py
import ufl
import adios4dolfinx
import pyvista as pv

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
        os.remove(os.getcwd() + "/Datastream.bp")
    except FileNotFoundError:
        pass
    ginput = read_geninput()

    filename = os.getcwd() + "/Datastream.xdmf"
    with io.XDMFFile(MPI.COMM_WORLD, filename, "a") as xdmf:
        xdmf.write_mesh(domain)
        V = fem.functionspace(domain, ("CG", 2, (1,)))
        for element in ginput["Material"]["Composition"].keys():
            u = fem.Function(V)
            u.name = "Composition/" + str(element)
            u.x.array[:] = np.ones(len(u.x.array[:]))*ginput["Material"]["Composition"][element]
            xdmf.write_function(u, 0.0)
        testu = fem.Function(V)
        testu.name = "test"
        testu.x.array[:] = 3.14159
        xdmf.write_function(testu, 0.0)

    print("Created datastream file")
def adjustdatastream(data, datapos="nodes", t_data=0.0):
    dataname = list(data.keys())[0]

    if not isinstance(data, dict):
        raise KeyError("datatype for datastream storage should be a dict")

    with io.XDMFFile(MPI.COMM_WORLD,"Datastream.xdmf", "r")as xdmf:
        domain = xdmf.read_mesh(name="Grid")

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
        #domain = adios4dolfinx.read_mesh("Datastream.xdmf", MPI.COMM_WORLD)
        #adios4dolfinx.read_function(function_filename.with_suffix(".bp"), u_in, time=0.0, name="Output")
        if dataname == "nodes":
            with io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf", "r") as xdmf:
                domain = xdmf.read_mesh(name="Grid")
            return domain.geometry.x[:]

        else:
            test = meshio.read("Datastream.xdmf")
            with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
                points, cells = reader.read_points_cells()
                if dataname == "nodes":
                    return points
                elif dataname == "elements":
                    return cells
                pd_list, cd_list, t_list = list(), list(), list()
                print("number of steps")
                print(reader.num_steps)
                for k in range(reader.num_steps):
                    print("Loop")
                    t, point_data, cell_data = reader.read_data(k)
                    print(point_data.keys())
                    if t == time and all_t == 0:
                        print(point_data[dataname])
                        return point_data[dataname]
                    t_list.append(t)
                    pd_list.append(point_data)
                    cd_list.append(cell_data)

            with io.XDMFFile(MPI.COMM_WORLD, "Datastream.bp", "r") as xdmf:
                domain = xdmf.read_mesh(name="Grid")
            V = fem.functionspace(domain, ("Lagrange", 2))
            u_in = fem.Function(V)

            print("Trying to read with Adios4")
            print("NOOOOOOOOOO")
            adios4dolfinx.read_function(os.getcwd() + "/Datastream.bp", u_in, time=0.0, name="Composition/Si")
            print(np.shape(u_in.x.array))
            return u_in.x.array[:]
    except:
        raise KeyError("!!!!!!!!")

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