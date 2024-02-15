import meshio
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
def readresultfile(filename, dataname, time=0):
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
        raise KeyError("Datastream "+str(dataname)+" doesn't exist in" + filename + "file. Data that exist is "
                       + str(pd_list[0].keys()))
def getaxisvalues_result(filename, dataname, time=0):
    node_y = readresultfile(filename, 'nodes')[:, 1]
    indx = np.where(node_y == 0)
    data = readresultfile(filename, dataname, time)[indx]
    x = readresultfile(filename, 'nodes')[:, 0][indx]

    # Sorting data in order from x=min(x) to x=max(x)
    indx = np.argsort(x)
    data = np.array(data)[indx]
    return data

def plotcompare(filenamelist, dataname, time=0):
    xvalues, yvalues = list(), list()
    for fn in filenamelist:
        xvalues.append(getaxisvalues_result(fn, 'nodes', time)[:, 0])
        yvalues.append(getaxisvalues_result(fn, dataname, time))
    for i in range(len(xvalues)):
        plt.plot(xvalues[i], yvalues[i],label=str(i))
    plt.legend()
    plt.show()
