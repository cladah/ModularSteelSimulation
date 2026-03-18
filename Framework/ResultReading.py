import meshio
import numpy as np
import json
def read_result_input(filename):
    f = open(filename,'r')
    data = json.load(f)
    return data
def getnames_results(filename):
    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        numsteps = reader.num_steps
        t, point_data, cell_data = reader.read_data(0)
        datanames = list(point_data.keys())
        if numsteps > 1:
            t, point_data, cell_data = reader.read_data(1)
            datanames = datanames + list(point_data.keys())
        return datanames
def get_times_results(filename, dataname):
    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        numsteps = reader.num_steps
        t, point_data, cell_data = reader.read_data(0)
        datanames = list(point_data.keys())
        if numsteps > 1:
            t, point_data, cell_data = reader.read_data(1)
            datanames = datanames + list(point_data.keys())
        return datanames
def read_results(filename, dataname, time=0):
    """
    Reading xdmf file and returning data for all points at a specific time step.

    :param filename: xdmf file where the results are stored
    :param dataname: The dataname of interest. Needs to be in xdmf keys
    :param time: time of interest. -1 is the last timestep. If nothing specified. First time step is used.
    :return:
    """

    try:
        pd_list, cd_list, t_list = list(), list(), list()
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            # Handling special data names
            if dataname == "nodes":
                return points
            elif dataname == "elements":
                return cells

            if time == -1: # Last timestep
                n = reader.num_steps-1
                t, point_data, cell_data = reader.read_data(n)
                if dataname not in point_data.keys():
                    available = list(point_data.keys())
                    raise KeyError(
                        f"Datastream '{dataname}' not found. Available datastreams: {available}")
                return point_data[dataname]


            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                t_list.append(t)
                pd_list.append(point_data)
                cd_list.append(cell_data)

                if t == time:
                    if dataname not in point_data:
                        available = list(point_data.keys())
                        raise KeyError(
                            f"Datastream '{dataname}' not found at time {t}. Available datastreams: {available}" )
                    return point_data[dataname]
            raise KeyError(f"Timestep {time} not found. Available timesteps: {t_list}")
    except meshio._exceptions.ReadError as e:
        available = list(pd_list[0].keys()) if pd_list else []
        raise KeyError(
            f"Error reading '{dataname}' from file '{filename}'. "
            f"Available datastreams: {available}"
        ) from e

def read_results_history(filename, dataname, xyz):
    """
    Reading xdmf file and returning data for a single point for all time steps.

    :param filename: xdmf file where the results are stored
    :param dataname: The dataname of interest. Needs to be in xdmf keys
    :param xyz: Cooridnate of interest
    :return t_list: Time array,
    :return pd_list: Data array
    """

    def find_nearest(array, value):
        ''' Find nearest value is an array '''
        idx = np.sqrt(np.power(np.subtract(array, value), 2).sum(1)).argmin()
        return idx

    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            indx = find_nearest(points, xyz)

            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                if k ==0.0:
                    continue
                else:
                    t, point_data, cell_data = reader.read_data(k)
                    t_list.append(t)
                    pd_list.append(point_data[dataname][indx])
            return t_list, pd_list
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))


def read_results_all(filename, points_xyz):
    """
    Returning all data from datastream file

    :param filename: xdmf file with data (str)
    :param points_xyz: Positions of points where data should be returned (list)
    :return: dict with np.arrays of data with datanames as keys
    """
    def find_nearest(array, value):
        ''' Find nearest value is an array '''
        idx = np.sqrt(np.power(np.subtract(array, value), 2).sum(1)).argmin()
        return idx

    try:
        pd_list, cd_list, t_list = list(), list(), list()
        pd_dict, cd_dict = dict(), dict()

        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            allpoints, cells = reader.read_points_cells()

            indxs = list()
            for point in points_xyz:
                indxs.append(find_nearest(allpoints, point))
            t, point_data0, cell_data = reader.read_data(0)
            if reader.num_steps > 1:
                t, point_data, cell_data = reader.read_data(1)

            for dataname in point_data0.keys():
                # Getting points on x-axis
                indx = np.where(allpoints[:,1] == 0)
                x = allpoints[:,0][indx]
                y = point_data0[dataname][indx]
                # Sorting datapoints
                indx = np.argsort(x)
                pd_dict[dataname] = [x[indx], y[indx]]
            if reader.num_steps == 1:
                return pd_dict
            for dataname in point_data.keys():
                pd_list = list()
                t_list = list()
                for k in range(reader.num_steps):
                    if k ==0.0:
                        continue
                    t, point_data, cell_data = reader.read_data(k)
                    t_list.append(t)
                    pd_list.append(point_data[dataname][indxs])
                pd_dict[dataname] = [t_list, pd_list]
            return pd_dict
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
def read_results_axis(filename, dataname, time=0):
    jsonfile = filename.replace(".xdmf", ".json")
    ginput = read_result_input(jsonfile)
    if ginput["Geometry"]["Type"] == "2D":
        node_y = read_results(filename, 'nodes')[:, 1]
        indx = np.where(node_y == 0)
        data = read_results(filename, dataname, time)[indx]
        x = read_results(filename, 'nodes')[:, 0][indx]

        # Sorting data in order from x=min(x) to x=max(x)
        indx = np.argsort(x)
        data = np.array(data)[indx]
        if dataname == "nodes":
            data = data[:, 0]
    elif ginput["Geometry"]["Type"] == "4PointBend":
        nodes = read_results(filename, 'nodes')
        if ginput["Geometry"]["dim"] == 3:
            #indx = np.where(np.isclose(nodes[:, 0], find_nearest(nodes[:, 0], 0.06)))
            indx = np.where(
                np.isclose(nodes[:, 0], 0.06) &
                np.isclose(nodes[:, 2], 0.0)
            )[0]
        else:
            print("Not implimented see geometry")
            raise KeyError()
        sort = np.argsort(nodes[indx][:,1])
        indx = np.take(np.array(indx), sort)
        data = read_results(filename, dataname, time)[indx]
        node_y = read_results(filename, 'nodes')[:, 1][indx]
        # Sorting data in order from x=min(x) to x=max(x)
        indx = np.argsort(node_y)
        data = np.array(data)[indx]
        if dataname == "nodes":
            data = data[:,1]
    else:
        raise KeyError("Geometry not implemented in read_results_axis")
    return data

def read_results_mesh(filename):
    import matplotlib as plt
    from mpi4py import MPI
    from dolfinx.io import gmshio, XDMFFile
    reader = meshio.xdmf.TimeSeriesReader(filename)
    num_steps = reader.num_steps
    points, cells = reader.read_points_cells()
    time_value, point_data, cell_data = reader.read_data(0)
    with XDMFFile(MPI.COMM_WORLD, filename, "r") as infile:
        # Read the mesh
        domain = infile.read_mesh()
    return domain