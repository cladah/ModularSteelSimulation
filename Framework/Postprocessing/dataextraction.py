from Framework.ResultReading import read_results_axis, getnames_results
import numpy as np
def DatastreamPlotting(dataname, t=0):
    import matplotlib.pyplot as plt
    filename = "Datastream.xdmf"
    print(getnames_results("Datastream.xdmf"))
    y = read_results_axis(filename, dataname, t)
    x = read_results_axis(filename, "nodes")
    plt.plot(x[:, 0], y)
    plt.show()

def ResultPlotting(filenames, dataname, point=[0., 0.], tid=-1):
    import matplotlib.pyplot as plt
    from Framework.HelpFile import get_plotlbls
    lbls = get_plotlbls(dataname)
    for filename in filenames:
        names = getnames_results(filename)
        print("Results in " + filename)
        print(names)
        tmpdata = read_results_axis(filename, dataname, tid)
        xyz = read_results_axis(filename, "nodes")
        plt.plot(xyz[:, 0]*1000, tmpdata)
    plt.legend(filenames)
    plt.title(dataname)
    plt.xlabel(lbls[0])
    plt.ylabel(lbls[1])
    plt.rcParams.update({'font.size': 30})
    plt.show()
    print("Done")

def export_data(filename, dataname, t=0):

    import pandas as pd
    y = read_results_axis(filename, dataname, t)
    x = read_results_axis(filename, "nodes")[:, 0]
    datadict = {"x": x, dataname: y}
    print(datadict)
    pd.DataFrame(datadict).to_csv("Resultfiles/TmpData.csv", index=False)

def xmdftesting():
    import meshio
    t_data = 0

    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()
        pd_list, cd_list, t_list = list(), list(), list()
        for k in range(reader.num_steps):
            t, point_data, cell_data = reader.read_data(k)
            t_list.append(t)
            pd_list.append(point_data)
            cd_list.append(cell_data)
    print(len(pd_list[0]["Composition/C"]))
    print(np.array(np.ones((len(pd_list[0]["Composition/C"]), 5))))

    data = {"ones":np.array(np.ones((len(pd_list[0]["Composition/C"]), 5)))}
    data = {"ones":np.array([[1,1,1,1,1,1] for i in range(len(pd_list[0]["Composition/C"]))])}

    if t_data not in t_list:
        t_list.append(t_data)
        pd_list.append(data)
        cd_list.append({})
    else:
        indx = t_list.index(t_data)
        for key in data.keys():
            pd_list[indx][key] = data[key]

    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(points, cells)
        for i in range(len(t_list)):
            writer.write_data(t=t_list[i], point_data=pd_list[i], cell_data=cd_list[i])