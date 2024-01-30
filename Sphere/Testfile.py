import meshio


def add_data_to_xdmf(filename, data, time):

    # Load the existing XDMF file if it exists
    try:
        mesh = meshio.read(filename)
    except FileNotFoundError:
        mesh = None

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

        # Add the new data for the specified time
        writer.write_data(time, point_data=data, cell_data=data)

    print(f"Data added for time {time} in {filename}")


def read_data_from_xdmf(filename, time):
    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            # Read mesh information
            points, cells, _, _ = reader.read_points_cells()

            # Read point data for the specified time
            data = reader.read_data(time)
            #data = reader.read_data(time=time)

            print(f"Data read for time {time} from {filename}")

            return points, cells, data

    except FileNotFoundError:
        print(f"File {filename} not found.")
        return None