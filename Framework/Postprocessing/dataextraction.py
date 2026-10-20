from Framework.ResultReading import read_results_axis, getnames_results, read_results, read_results_all, read_line_data
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
        plt.plot(xyz[:]*1000, tmpdata)
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
    x = read_results_axis(filename, "nodes")

    datadict = {"x": x}
    for dataname in datanames:
        print(dataname)
        y = read_results_axis(filename, dataname, t)
        if len(np.shape(y)) == 1:
            datadict[dataname] = y
    datapd = pd.DataFrame(datadict).sort_values('x')

    num_rows = len(datapd)
    """
    # Linearly increasing
    roll = 50
    datapd.rolling(roll, 1)
    indices = np.geomspace(1, num_rows - 1, 50).astype(int)
    indices = np.append([0], indices)
    indices = np.unique(indices)
    datapd = datapd.iloc[indices, :]
    print(datapd)
    print(np.shape(datapd))
    """
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
    print(np.shape(data))
    #if np.shape(data)[1] == 6:
    #    plt.legend(["x","y","z","xy","yz","xz"])
    plt.show()

def plotting_over_line(filename="Datastream.xdmf"):
    available_fields = getnames_results(filename)
    for i, f in enumerate(available_fields): print(f" [{i}] {f}")
    field_idx = int(input("Choice: "))
    selected_field = available_fields[field_idx]
    for dist in np.linspace(0.0,0.002,10):
        x, data = read_line_data(filename, selected_field, [0.0, dist, 0.0],[1.0, 0.0, 0.0], mirror_x=0.06)
        plt.plot(x,data[:, 0]/1e6,label=f"Carbon content {np.round(dist*1000,2)} mm")
    plt.xlabel("Distance along beam (mm)")
    plt.ylabel("Stress along beam (MPa)")
    plt.legend()
    plt.show()





def _select_component(data: np.ndarray, field_name: str) -> np.ndarray:
    STRESS_COMPONENTS = ("x", "y", "z", "xz", "yz", "xy")
    """Return a scalar array for PyVista coloring."""
    data = np.asarray(data)

    if data.ndim == 1:
        return data

    if data.ndim < 2:
        raise ValueError(f"Unexpected data shape: {data.shape}")

    # Flatten tensor dimensions, preserving the point dimension.
    components = data.reshape(data.shape[0], -1)
    n_components = components.shape[1]

    print(f"\nField '{field_name}' has shape {data.shape}.")
    print("Select component:")
    print(" [6] Magnitude / Frobenius norm")

    if n_components == 3:
        labels = ("x", "y", "z")
    elif n_components == 6 and field_name in {
        "Stress", "Strain_el"
    }:
        labels = STRESS_COMPONENTS
    else:
        labels = tuple(str(i) for i in range(n_components))

    for i, label in enumerate(labels):
        print(f" [{i}] {label}")

    while True:
        try:
            choice = int(input("Component: "))
            if choice == 6:
                return np.linalg.norm(components, axis=1)
            if 0 <= choice < n_components:
                return components[:, choice]
        except ValueError:
            pass

        print(f"Enter 0–{n_components - 1}, or 6 for magnitude.")
def _build_grid(
    points: np.ndarray,
    cells: list,
) -> pv.UnstructuredGrid:
    """Convert all supported meshio cell blocks to a PyVista grid."""
    CELL_TYPE_MAP = {
        "vertex": 1,
        "line": 3,
        "triangle": 5,
        "quad": 9,
        "tetra": 10,
        "hexahedron": 12,
        "wedge": 13,
        "pyramid": 14,
        "triangle6": 22,
        "quad8": 23,
        "tetra10": 24,
        "hexahedron20": 25,
    }

    points = np.asarray(points, dtype=float)

    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ValueError(
            f"Expected 2D or 3D points; got shape {points.shape}."
        )

    if points.shape[1] == 2:
        points = np.column_stack(
            (points, np.zeros(points.shape[0]))
        )

    connectivity_blocks = []
    cell_type_blocks = []

    for block in cells:
        vtk_type = CELL_TYPE_MAP.get(block.type)

        if vtk_type is None:
            print(f"Skipping unsupported cell type: {block.type}")
            continue

        connectivity = np.asarray(block.data, dtype=np.int64)

        if connectivity.ndim != 2:
            raise ValueError(
                f"Invalid connectivity for cell type {block.type}."
            )

        if connectivity.size and (
            connectivity.min() < 0
            or connectivity.max() >= len(points)
        ):
            raise ValueError(
                f"Invalid point indices in {block.type} cells."
            )

        n_points_per_cell = connectivity.shape[1]

        vtk_cells = np.column_stack(
            (
                np.full(len(connectivity), n_points_per_cell),
                connectivity,
            )
        )

        connectivity_blocks.append(vtk_cells.ravel())
        cell_type_blocks.append(
            np.full(len(connectivity), vtk_type, dtype=np.uint8)
        )

    if not connectivity_blocks:
        raise ValueError("No supported cells were found in the mesh.")

    # Each block can have a different number of nodes per cell.
    vtk_connectivity = np.concatenate(connectivity_blocks)
    vtk_cell_types = np.concatenate(cell_type_blocks)

    return pv.UnstructuredGrid(
        vtk_connectivity,
        vtk_cell_types,
        points,
    )
def interactive_xdmf_plotter(
    filename: str = "Datastream.xdmf",
    *,
    warp: bool = False,
    warp_factor: float = 2.0,
    show_edges: bool = False,
    show_contours: bool = False,
    contour_count: int = 15,
    window_size: tuple[int, int] = (1600, 1000),
) -> None:
    """
    Interactively select an XDMF time step and point-data field,
    then visualize the result using PyVista.

    Parameters
    ----------
    filename:
        Path to the meshio-compatible XDMF time-series file.
    warp:
        Deform the mesh using the full Displacement vector, if present.
    warp_factor:
        Scale factor for displacement-based deformation.
    show_edges:
        Display element edges.
    show_contours:
        Overlay scalar contours on the external surface.
    contour_count:
        Number of contour levels.
    window_size:
        PyVista render-window dimensions.
    """
    print(f"\n--- Opening {filename} ---")

    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            num_steps = reader.num_steps

            if num_steps < 1:
                print("No time steps found.")
                return

            # Inspect the first time step for available fields.
            _, first_point_data, _ = reader.read_data(0)
            initial_fields = list(first_point_data)

            if not initial_fields:
                print("No point-data fields found.")
                return

            print(f"\nTime steps: {num_steps}")
            print("Fields in step 0:", ", ".join(initial_fields))

            # Select time step.
            while True:
                try:
                    step_idx = int(
                        input(f"Select time step (0–{num_steps - 1}): ")
                    )
                    if 0 <= step_idx < num_steps:
                        break
                except ValueError:
                    pass

                print("Invalid time-step index.")

            time_value, point_data, _cell_data = reader.read_data(step_idx)

            if not point_data:
                print(f"No point data found at time step {step_idx}.")
                return

            fields = list(point_data)

            print("\nAvailable fields:")
            for i, name in enumerate(fields):
                print(f" [{i}] {name}")

            # Select field.
            while True:
                try:
                    field_idx = int(input("Select field: "))
                    if 0 <= field_idx < len(fields):
                        break
                except ValueError:
                    pass

                print("Invalid field index.")

            selected_field = fields[field_idx]
            raw_data = np.asarray(point_data[selected_field])

            if raw_data.shape[0] != len(points):
                raise ValueError(
                    f"Field '{selected_field}' has {raw_data.shape[0]} "
                    f"entries, but the mesh has {len(points)} points. "
                    "This function expects point data."
                )

            scalar_data = _select_component(raw_data, selected_field)

            if not np.all(np.isfinite(scalar_data)):
                print(
                    "Warning: the selected scalar field contains "
                    "NaN or infinite values."
                )

            print(
                f"\nLoading '{selected_field}' at "
                f"t = {float(time_value):.6g}"
            )

            # Preserve the original displacement vector independently
            # of the selected field and its scalar component.
            displacement = point_data.get("Displacement")

            # Construct the grid while the mesh data is available.
            grid = _build_grid(points, cells)

    except (OSError, ValueError, KeyError) as exc:
        print(f"Could not load XDMF data: {exc}")
        return

    # Store a scalar array for coloring.
    scalar_name = f"{selected_field}"
    grid.point_data[scalar_name] = scalar_data
    if scalar_name == "Temperature":
        grid.point_data[scalar_name] -= 273.15

    # Warp using the complete displacement vector, not a scalar component.
    if warp:
        if displacement is None:
            print(
                "Warning: Displacement is unavailable; "
                "plotting the undeformed mesh."
            )
        else:
            displacement = np.asarray(displacement)

            if (
                displacement.ndim != 2
                or displacement.shape[0] != len(points)
                or displacement.shape[1] not in (2, 3)
            ):
                print(
                    "Warning: Displacement must have shape (N, 2) "
                    "or (N, 3); skipping deformation."
                )
            else:
                if displacement.shape[1] == 2:
                    displacement = np.column_stack(
                        (
                            displacement,
                            np.zeros(len(displacement)),
                        )
                    )

                # Avoid overwriting the selected scalar field.
                grid.point_data["_warp_displacement"] = displacement
                grid = grid.warp_by_vector(
                    "_warp_displacement",
                    factor=warp_factor,
                )

    # Select a colormap based on the field name.
    name_lower = selected_field.lower()
    if "stress" in name_lower:
        cmap = "turbo"
    elif "temp" in name_lower:
        cmap = "inferno"
    elif "displacement" in name_lower:
        cmap = "viridis"
    else:
        cmap = "viridis"

    # Plot.
    plotter = pv.Plotter(window_size=window_size)

    plotter.add_text(
        f"Field: {selected_field} | "
        f"Time: {float(time_value):.4g}",
        font_size=12,
    )
    grid = grid.scale(1000)
    plotter.add_mesh(
        grid,
        scalars=scalar_name,
        cmap=cmap,
        show_edges=show_edges,
        scalar_bar_args={"title": selected_field},
    )

    if show_contours:
        try:
            surface = grid.extract_surface().triangulate()
            values = np.asarray(surface.point_data[scalar_name])
            lo, hi = np.nanmin(values), np.nanmax(values)

            if np.isfinite(lo) and np.isfinite(hi) and lo < hi:
                contours = surface.contour(
                    isosurfaces=np.linspace(
                        lo, hi, contour_count + 2
                    )[1:-1],
                    scalars=scalar_name,
                )
                plotter.add_mesh(
                    contours,
                    color="white",
                    line_width=1.5,
                )
            else:
                print("Skipping contours: scalar range is constant.")

        except (ValueError, RuntimeError) as exc:
            print(f"Contour generation skipped: {exc}")

    plotter.enable_parallel_projection()
    plotter.camera_position = "iso"
    plotter.camera.up = (0, 1, 0)
    plotter.reset_camera()
    plotter.show_grid(
        xtitle="X [mm]",
        ytitle="Y [mm]",
        ztitle="Z [mm]",
        font_size=10,
        fmt="%.0f"
    )
    plotter.show_axes()
    plotter.show()
def interactive_xdmf_plotter_old(filename="Datastream.xdmf"):
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
            t, point_data, cell_data = reader.read_data(step_idx)
            available_fields = list(point_data.keys())
            print("\nSelect Field:")
            for i, f in enumerate(available_fields): print(f" [{i}] {f}")
            field_idx = int(input("Choice: "))
            selected_field = available_fields[field_idx]
            data_array = point_data[selected_field]
        except (ValueError, IndexError):
            print("Invalid selection. Exiting.")
            return

        # 3. Read the specific data

        #warp_array = point_data["Displacement"]
        print(f">>> Loading {selected_field} at t = {t:.4f}s")
    selected_axis = 0
    if selected_field in ["Stress", "Strain_el"]:
        print("Select axis:"
              "\n[6] Norm"
              "\n[0] x"
              "\n[1] y"
              "\n[2] z"
              "\n[3] xz"
              "\n[4] yz"
              "\n[5] xy")
        selected_axis = input("Choice: ")
        if not selected_axis.isdigit():
            print("Invalid input. Please enter a number.")
            selected_axis = int(input("Choice: "))
        selected_axis = int(selected_axis)
        if selected_axis in [0, 1, 2, 3, 4, 5]:
            print(np.shape(data_array))
            data_array = data_array[:, selected_axis]
    elif selected_field in ["Displacement"]:
        print(data_array)
        print(np.shape(data_array))
        print("Select axis:"
              "\n[6] Norm"
              "\n[0] x"
              "\n[1] y"
              "\n[2] z")
        selected_axis = int(input("Choice: "))
        if selected_axis in [0, 1, 2]:
            print(np.shape(data_array))
            data_array = data_array[:, selected_axis]

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
        #grid.point_data["Displacement"] = warp_array

    # 6. Plotting
    p = pv.Plotter(window_size=[2560,1440])
    p.add_text(f"Field: {selected_field} | Time: {t:.4f}s", font_size=6)

    # Auto-color logic
    cmap = "turbo" if "Stress" in selected_field else "inferno" if "Temp" in selected_field else "viridis"
    """
    axis_data = read_results_axis(filename, selected_field)
    if len(np.shape(axis_data)) > 1:
        if selected_axis != 6:
            axis_data = axis_data[:, selected_axis]
        else:
            axis_data = np.abs(axis_data)
    vmax = np.max(axis_data)
    vmin = np.min(axis_data)
    """
    displacement_key = "Displacement"
    warp_factor = 2.0
    print(grid.point_data)
    if displacement_key in grid.point_data:
        print("Deforming")
        # Warping before reflection ensures symmetry is perfectly preserved
        #grid = grid.warp_by_vector(displacement_key, factor=warp_factor)
    else:
        print(f">>> Warning: '{displacement_key}' not found in grid point data. Skipping warp.")

    if False:
        reflected_grid = grid.reflect(normal=(1, 0, 0), point=(0.06, 0, 0))
        grid = grid + reflected_grid
        reflected_grid = grid.reflect(normal=(0, 0, 1), point=(0, 0, 0))
        grid = grid + reflected_grid



    p.add_mesh(grid, scalars=selected_field, cmap=cmap, show_edges=False)

    try:
        external_surface = grid.extract_surface()
        contours = external_surface.triangulate().contour(
            isosurfaces=15,
            scalars=selected_field
        )

        # 4. Add the lines to the plotter
        #p.add_mesh(contours, color="white", line_width=1, label="Boundary Contours")

    except Exception as e:
        print(f">>> Contour generation skipped: {e}")
    p.enable_parallel_projection()
    p.camera.up = (0, 1, 0)

    # Position the camera in the (-x, +y, -z) quadrant
    # This makes +x point down-right and +z point up-right
    p.camera.position = (-3, 1, -4)

    # Point the camera at the center of your data
    p.camera.focal_point = (0, 0, 0)
    p.reset_camera()
    #p.view_xy()
    #p.show_grid(font_size=10)

    print(">>> Launching PyVista window...")
    p.save_graphic("Resultfiles/pyvista.pdf")
    p.show(screenshot="Resultfiles/pyvista.png")



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