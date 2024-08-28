import meshio
import numpy as np
def getnames_results(filename):
    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            datanames = point_data.keys()
            return datanames

def read_results(filename, dataname, time=0):
    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
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
            if time == -1:
                maxindx = t_list.index(np.max(t_list))
                return pd_list[maxindx][dataname]
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))
    return

def read_results_history(filename, dataname, xyz):
    def find_nearest(array, value):
        ''' Find nearest value is an array '''
        idx = (np.abs(array - value)).argmin()
        return idx

    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            points()

            pd_list, cd_list, t_list = list(), list(), list()
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                t_list.append(t)
                pd_list.append(point_data[dataname])
            return t_list, pd_list
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))
    return