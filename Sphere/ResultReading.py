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
    return
def read_results_all(filename, xyz):
    def find_nearest(array, value):
        ''' Find nearest value is an array '''
        idx = np.sqrt(np.power(np.subtract(array, value), 2).sum(1)).argmin()
        return idx

    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            indxs = list()
            for point in xyz:
                indxs.append(find_nearest(points, point))
            t, point_data0, cell_data = reader.read_data(0)
            t, point_data, cell_data = reader.read_data(1)
            pd_list, cd_list, t_list = list(), list(), list()
            pd_dict, cd_dict = dict(), dict()
            for dataname in point_data0.keys():
                # Getting points on x-axis
                indx = np.where(points[:,1] == 0)
                x = points[:,0][indx]
                y = point_data0[dataname][indx]
                # Sorting datapoints
                indx = np.argsort(x)
                pd_dict[dataname] = [x[indx], y[indx]]
            for dataname in point_data.keys():
                pd_list = list()
                t_list = list()
                for k in range(reader.num_steps):
                    if k ==0.0:
                        continue
                    else:
                        t, point_data, cell_data = reader.read_data(k)
                        t_list.append(t)
                        pd_list.append(point_data[dataname][indxs])
                pd_dict[dataname] = [t_list, pd_list]
            return pd_dict
    except meshio._exceptions.ReadError:
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in datastream file. Data that exist is "
                       + str(pd_list[0].keys()))
    return