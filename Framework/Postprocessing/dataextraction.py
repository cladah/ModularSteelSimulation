from Framework.ResultReading import read_results_axis, getnames_results, read_results, read_results_all
import numpy as np
import dolfinx
from dolfinx.io import XDMFFile
from mpi4py import MPI
import pyvista as pv
import meshio
import matplotlib.pyplot as plt

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

def export_data(filename, datanames, t=0, n=1):
    print("\n---------------------------------------------------------------------\n")
    print("Exporting data!")
    if datanames == "All":
       datanames = getnames_results(filename)
    elif type(datanames) != list:
        datanames = [datanames]

    import pandas as pd
    x = read_results(filename, "nodes")[:, 0]
    #x = read_results_axis(filename, "nodes")[:, 0]
    datadict = {"x": x}
    for dataname in datanames:
        print(dataname)
        y = read_results(filename, dataname, t)
        #y = read_results_axis(filename, dataname, t)
        if len(np.shape(y)) == 1:
            datadict[dataname] = y
    datapd = pd.DataFrame(datadict).sort_values('x')

    num_rows = len(datapd)
    # Linearly increasing
    roll = 50
    datapd.rolling(roll, 1)
    indices = np.geomspace(1, num_rows - 1, 50).astype(int)
    indices = np.append([0], indices)
    indices = np.unique(indices)
    datapd = datapd.iloc[indices, :]
    datapd.to_csv("Resultfiles/TmpData.csv", index=False)
    print("\n---------------------------------------------------------------------\n")
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

def plotting_over_axis(filename="Datastream.xdmf"):
    available_fields = getnames_results(filename)
    for i, f in enumerate(available_fields): print(f" [{i}] {f}")
    field_idx = int(input("Choice: "))
    selected_field = available_fields[field_idx]
    data = read_results_axis(filename, selected_field)
    x = read_results_axis(filename, "nodes")

    plt.plot(x,data)
    plt.show()


def interactive_datastream_plotter(filename="Datastream.xdmf"):
    """
    Uses meshio TimeSeriesReader to navigate through time steps
    and plot specific functions via PyVista.
    """
    print(f"\n--- Opening {filename} with TimeSeriesReader ---")

    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        # 1. Get basic mesh info
        points, cells = reader.read_points_cells()
        num_steps = reader.num_steps

        # Read the first step just to see what data fields exist
        _, point_data, _ = reader.read_data(0)
        available_fields = list(point_data.keys())

        if not available_fields:
            print(">>> No data fields found in the file.")
            return

        # 2. User Menu for Step and Field
        print(f"\nFile contains {num_steps} time steps.")
        print("Available Fields:", ", ".join(available_fields))

        try:
            step_idx = int(input(f"Select Time Step (0 to {num_steps - 1}): "))
            print("\nSelect Field:")
            for i, f in enumerate(available_fields): print(f" [{i}] {f}")
            field_idx = int(input("Choice: "))
            selected_field = available_fields[field_idx]
        except (ValueError, IndexError):
            print("Invalid selection. Exiting.")
            return

        # 3. Read the specific data
        t, point_data, cell_data = reader.read_data(step_idx)
        data_array = point_data[selected_field]
        print(f">>> Loading {selected_field} at t = {t:.4f}s")

    # 4. Construct PyVista Grid
    cell_map = {
        "triangle": 5,  # Linear Triangle (3-node)
        "triangle6": 22,  # Quadratic Triangle (6-node)
        "quad": 9,  # Linear Quad (4-node)
        "quad8": 23,  # Quadratic Quad (8-node)
        "tetra": 10,  # Linear Tetra (4-node)
        "tetra10": 24,  # Quadratic Tetra (10-node)
        "hexahedron": 12,  # Linear Hex (8-node)
        "hexahedron20": 25  # Quadratic Hex (20-node)
    }
    cell_conn = None
    vtk_type = None

    for block in cells:
        if block.type in cell_map:
            cell_conn = block.data
            vtk_type = cell_map[block.type]
            break
    if cell_conn is None:
        print(f"Unsupported cell type in mesh: {[b.type for b in cells]}")
        return

    # Format points for 3D engine
    if points.shape[1] == 2:
        points = np.hstack([points, np.zeros((points.shape[0], 1))])

    # Build the connectivity array: [n_pts, p1, p2, p3, ...]
    cells_pv = np.hstack([np.full((len(cell_conn), 1), cell_conn.shape[1]), cell_conn])
    grid = pv.UnstructuredGrid(cells_pv, np.full(len(cell_conn), vtk_type), points)

    # 5. Handle Vector Magnitude (if needed)
    if data_array.ndim > 1:
        # Calculate magnitude for coloring
        mag = np.linalg.norm(data_array, axis=1)
        grid.point_data[selected_field] = mag
        # Store actual vectors for potential warping/glyphs
        grid.point_data[f"{selected_field}_vec"] = np.hstack([data_array, np.zeros((data_array.shape[0], 1))]) if \
        data_array.shape[1] == 2 else data_array
    else:
        grid.point_data[selected_field] = data_array

    # 6. Plotting
    p = pv.Plotter()
    p.add_text(f"Field: {selected_field} | Time: {t:.4f}s", font_size=6)

    # Auto-color logic
    cmap = "turbo" if "Stress" in selected_field else "inferno" if "Temp" in selected_field else "viridis"

    p.add_mesh(grid, scalars=selected_field, cmap=cmap, show_edges=True)
    #p.add_scalar_bar()
    p.view_xy()
    p.show_grid(font_size=10)

    print(">>> Launching PyVista window...")
    p.show()



def testread(filename):
    t = 0.0
    # 1. Get the list of available data names
    available_names = getnames_results(filename)

    if not available_names:
        print("No data names found in the file.")
        return

    # 2. Display the menu
    print("\n--- Available Data Names ---")
    for i, name in enumerate(available_names, 1):
        print(f"{i}. {name}")
    print("----------------------------")

    # 3. Get user input
    choice = input("\nEnter the name (or number) you want to print: ").strip()

    # 4. Handle numerical selection or direct name entry
    if choice.isdigit():
        idx = int(choice) - 1
        if 0 <= idx < len(available_names):
            dataname = available_names[idx]
        else:
            print("Invalid selection.")
            return
    else:
        dataname = choice

    # 5. Read and print the results
    if dataname in available_names:
        print(f"\nResults for '{dataname}' at t={t}:")
        data = read_results_axis(filename, dataname, t)
        nodes = read_results_axis(filename, "nodes")
        import matplotlib.pyplot as plt
        print(nodes)
        print(data)
        plt.plot(nodes, data)
        plt.show()
        print(data)
    else:
        print(f"Error: '{dataname}' not found in results.")
    return
    import h5py
    import pyvista as pv
    xdmf = pv.read(filename)
    if isinstance(xdmf, pv.MultiBlock):
        print("Multiblock")
        print(xdmf.keys())
        timestep = xdmf[1]
    else:
        timestep = xdmf
    print(timestep.GetPointData())
    for method in dir(timestep):
        if callable(getattr(timestep, method)):
            #print(method)
            pass
    #print([method_name for method_name in dir(timestep)
    # if callable(getattr(timestep, method_name))])

    return timestep[dataname]
    if dataname in mesh.point_data:
        data = mesh.point_data[dataname]
    elif dataname in mesh.cell_data:
        data = mesh.cell_data[dataname]
    else:
        raise KeyError(f"Array '{dataname}' not found. "
                       f"Available arrays: {list(mesh.point_data.keys()) + list(mesh.cell_data.keys())}")

    return data


    with h5py.File(filename, "r") as f:
        for test in f:
            print(test)
        if dataname not in f:
            raise KeyError(f"Dataset '{filename}' not found in {filename}. "
                           f"Available keys: {list(f.keys())}")
        data = f[dataname][()]  # Read dataset into numpy array
    print(data)
    return data