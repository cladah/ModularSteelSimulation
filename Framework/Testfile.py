import meshio
from Datastream_file import *
from HelpFile import read_input
from Modulefiles.TTTdiagram_file import runTTTcalc
def add_data_to_xdmf(filename, data, time):

    # Load the existing XDMF file if it exists
    try:
        mesh = meshio.read(filename)
    except FileNotFoundError:
        mesh = None

    meshio.write("Resultfiles/Test.xdmf", mesh)
    # Create or update the TimeSeriesWriter
    with meshio.xdmf.TimeSeriesWriter("Resultfiles/Test.xdmf") as writer:
        if mesh is not None:
            # Add existing mesh data
            writer.write_points_cells(mesh.points, mesh.cells)

            # Add existing point data
            for name, array in mesh.point_data.items():
                writer.write_data(name, array)

            # Add existing field data
            for name, array in mesh.field_data.items():
                writer.write_data(name, array)
        if data:
            # Add the new data for the specified time
            writer.write_data(time, point_data=data, cell_data=data)

    print(f"Data added for time {time} in {filename}")


def read_data_from_xdmf(filename, time):
    pass

def testdatastream():
    print(gethistoryvalues("T", 0))

def createTTTdiagram_loop():
    data = read_input()
    composition = data["Material"]["Composition"]
    compositions = list()
    Cvariation = np.arange(0, 2.1, 0.1)
    Nvariation = np.arange(0, 1.6, 0.1)
    for i in Cvariation:
        for j in Nvariation:
            tmpcomp = composition.copy()
            tmpcomp["C"] = i
            tmpcomp["N"] = j
            compositions.append(tmpcomp)
    for comp in compositions:
        runTTTcalc(comp)
