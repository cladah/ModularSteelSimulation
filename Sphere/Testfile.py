import meshio

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