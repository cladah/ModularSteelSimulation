
import numpy as np
import meshio
from HelpFile import read_input, createinputcache, createresultinput, read_geninput, read_modinput
import json
import os
import pathlib
from pathlib import Path
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

def createdatastreaminput():
    ginput = read_geninput()
    f = open("Datastream.json", "w")
    totinput = dict()
    totinput.update(ginput)
    for file in ginput["Inputs"]:
        minput = read_modinput("Inputs/" + ginput["InputDirectory"] + "/" + file + ".json")
        totinput.update(minput)
    json.dump(totinput, f, indent=2)
    f.close()


def createdatastream(domain, name="Datastream"):
    """
    Initializes XDMF, H5, and JSON files for a simulation domain.
    """
    # 1. Setup paths using pathlib
    base_path = Path.cwd()
    files = {ext: base_path / f"{name}.{ext}" for ext in ["xdmf", "h5", "json"]}

    # 2. Cleanup existing files
    for path in files.values():
        if path.exists():
            path.unlink()
            print(f"Removed existing file: {path.name}")

    # 3. Process Input Data
    try:
        ginput = read_geninput()
    except Exception as e:
        raise RuntimeError(f"Failed to read general input: {e}")

    # Merge inputs efficiently
    total_input = ginput.copy()
    input_dir = Path("Inputs") / ginput.get("InputDirectory", "")

    for file_name in ginput.get("Inputs", []):
        mod_path = input_dir / f"{file_name}.json"
        try:
            mod_input = read_modinput(str(mod_path))
            total_input.update(mod_input)
        except Exception as e:
            print(f"Warning: Could not read module input {mod_path}: {e}")

    # Save merged JSON
    with open(files["json"], "w") as f:
        json.dump(total_input, f, indent=2)

    # 4. Initialize XDMF Mesh
    with meshio.xdmf.TimeSeriesWriter(files["xdmf"]) as writer:
        cell_types = list(domain.cells_dict.keys())
        if not cell_types:
            raise ValueError("Domain contains no cells.")

        # Flexibility: Use all coordinates (3D) or slice for 2D
        # Most solvers expect 3D points even for 2D meshes (z=0)
        points = domain.points

        # Handle mixed meshes or just take the primary type
        main_cell_type = cell_types[0]
        cells = {main_cell_type: domain.get_cells_type(main_cell_type)}

        writer.write_points_cells(points, cells=cells)

        # Initialize the time series at t=0
        writer.write_data(0.0, point_data={}, cell_data={})

    # 5. Initialize Composition Data
    composition = ginput.get("Material", {}).get("Composition", {})
    for element, value in composition.items():
        # Create an array of constant values across the points
        initial_values = np.full(len(domain.points), value)
        adjustdatastream({f"Composition_{element}": initial_values}, datapos="nodes")

    print(f"Successfully created datastream: {name}")

def adjustdatastream(data, filename="Datastream.xdmf", datapos="nodes", t_data=0.0):
    """
    Updates or appends data to an XDMF time series.

    Args:
        data (dict): Dictionary of numpy arrays to store.
        filename (str): Path to the .xdmf file.
        datapos (str): "nodes" for point_data, "cells" for cell_data.
        t_data (float): The timestamp to target.
    """
    if not isinstance(data, dict):
        raise TypeError(f"Data must be a dict, got {type(data).__name__}")

    file_path = Path(filename)
    if not file_path.exists():
        raise FileNotFoundError(f"Could not find {filename}. Create the mesh first.")

    # 1. Read existing state
    all_steps = []
    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        n_steps = reader.num_steps

        for k in range(n_steps):
            t, p_data, c_data = reader.read_data(k)
            all_steps.append({"t": t, "point_data": p_data, "cell_data": c_data})

    # 2. Update or Append logic
    found_time = False
    for step in all_steps:
        if np.isclose(step["t"], t_data):
            # Update existing timestep
            target_key = "point_data" if datapos == "nodes" else "cell_data"
            step[target_key].update(data)
            found_time = True
            break

    if not found_time:
        # Create a new timestep entry
        new_step = {
            "t": t_data,
            "point_data": data if datapos == "nodes" else {},
            "cell_data": data if datapos == "cells" else {}
        }
        all_steps.append(new_step)
        # Ensure steps stay sorted by time
        all_steps.sort(key=lambda x: x["t"])

    # 3. Write everything back
    with meshio.xdmf.TimeSeriesWriter(filename) as writer:
        writer.write_points_cells(points, cells)
        for step in all_steps:
            writer.write_data(
                step["t"],
                point_data=step["point_data"],
                cell_data=step["cell_data"]
            )

    action = "Updated" if found_time else "Appended new"
    print(f"{action} timestep {t_data} in {filename}")
def getnamesdatastream()-> list[str]:
    data_names = []
    xdmf_path = os.path.join(os.getcwd(), "Datastream.xdmf")
    h5_path = os.path.join(os.getcwd(), "Datastream.h5")

    # Checking if there is a datastream
    if not os.path.exists("Datastream.xdmf") or not os.path.exists(h5_path):
        print(f"Error: One or both files not found ({xdmf_path}, {h5_path}).")
        return data_names

    with io.XDMFFile(MPI.COMM_WORLD, "Datastream.xdmf") as xdmf:
        mesh = xdmf.read_mesh()

        with h5py.File(h5_path, "r") as f:
            if "Function" in f:
                function_group = f["Function"]
                for key in function_group.keys():
                    data_names.append(key)
                    print(f"  Found data field: **{key}**")
            else:
                print(f"Warning: 'Function' group not found in {h5_path}.")
    return data_names

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
                #print(f"  Step {k}: time={t}, point_data={list(point_data.keys())}")
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
