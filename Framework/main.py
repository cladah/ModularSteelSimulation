import time
import tkinter
import tkinter as tk
import sys
import json
from tkinter import filedialog
from pathlib import Path
from CGUImodule import MainApp
from Result_GUI import Result_MainApp, Compare_MainApp
from Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import createinputcache, read_geninput, reset_output, get_file_path, checkDatabase
import customtkinter as ctk
from Postprocessing.dataextraction import ResultPlotting, export_data, interactive_datastream_plotter, plotting_over_axis
import MSSModule as MSS

def progressmonitor(tid, module):
    """
    :param tid:
    :param module:
    """
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)

def GUI():
    """
    Using Tkinter to run the GUI
    """
    print("Opening GUI window...")
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()

def Result_GUI_comp(filenames, dataname):
    """
            Using Tkinter to run the GUI
        """
    print("Opening result window...")
    ctk.set_appearance_mode("dark")
    app = Compare_MainApp(filenames, dataname)
    app.mainloop()

def Result_GUI_show(filename):
    """
        Using Tkinter to run the GUI
    """
    if filename == "":
        root = tkinter.Tk()
        root.withdraw()
        filename = tkinter.filedialog.askopenfilename(filetypes=(("xdmf files", "*.xdmf"),))
    print("Opening result window...")
    ctk.set_appearance_mode("dark")
    app = Result_MainApp(filename)
    app.mainloop()

def vtxfile_test():
    from adios2 import FileReader

    with FileReader("FeniCSx/beam_stress.bp") as s:
        # inspect variables
        vars = s.available_variables()
        for name, info in vars.items():
            print("variable_name: " + name, end=" ")
            for key, value in info.items():
                print("\t" + key + ": " + value, end=" ")
            print()
        print()
    pass
def Test():
    import gmsh
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    import numpy as np

    # --- 1. Gmsh Setup & Mesh Generation ---
    gmsh.initialize()
    gmsh.model.add("QuarterCirclePlot")

    r = 1.0
    n_nodes = 25
    progression = 0.8  # Adjust > 1 for density at the edge

    # Geometry
    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint(r, 0, 0)
    p3 = gmsh.model.geo.addPoint(0, r, 0)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addCircleArc(p2, p1, p3)
    l3 = gmsh.model.geo.addLine(p3, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3])
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Transfinite constraints (density away from p1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, n_nodes, "Power", progression)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, n_nodes, "Power", 1 / progression)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, n_nodes)
    gmsh.model.geo.mesh.setTransfiniteSurface(s)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # --- 2. Extract Data for Plotting ---
    # Get node coordinates
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = coords.reshape((-1, 3))
    x = nodes[:, 0]
    y = nodes[:, 1]

    # Get element connectivity (Triangles)
    _, element_tags, node_connectivity = gmsh.model.mesh.getElements(2)
    # Gmsh returns 1D array of node tags; for triangles, we group by 3
    # Note: Tags are 1-based, Python is 0-based, so subtract 1
    triangles = (node_connectivity[0].reshape((-1, 3)) - 1)

    # --- 3. Matplotlib Plotting ---
    plt.figure(figsize=(8, 8))
    plt.gca().set_aspect('equal')

    # Create the triangulation object
    triangulation = tri.Triangulation(x, y, triangles)

    # Plot the mesh edges
    plt.triplot(triangulation, 'b-', lw=0.5)

    # Optional: Add some styling
    plt.title(f"Transfinite Mesh (Progression: {progression})", fontsize=14)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True, linestyle='--', alpha=0.6)

    print("Mesh generated. Close the plot window to finalize Gmsh.")
    plt.show()

    # Cleanup
    gmsh.finalize()

def main():
    menu = (
        "\n--- Simulation Manager ---\n"
        "1 - Run Simulation\n"
        "2 - Run using GUI\n"
        "3 - View Results\n"
        "4 - Export Data\n"
        "5 - Compare Files\n"
        "6 - Debug/Test Read\n"
        "7 - Plot value over axis\n"
        "8 - Check CHALPHAD database\n"
        "9 - Optimization\n"
        "10 - Testing\n"
        "0 - Exit\n"
        "Selection: "
    )
    while True:
        user_input = input(menu).strip()
        print(" ")
        if not user_input.isdigit():
            print("Invalid input. Please enter a number.")
            continue

        choice = int(user_input)
        match choice:
            case 1:
                done = MSS.run_simulation()
                if done:
                    interactive_datastream_plotter()
            case 2:
                GUI()
            case 3:
                Result_GUI_show("")
            case 4:
                filename = get_file_path("Select XDMF to Export")
                if filename:
                    export_data(filename, ["Composition/C", "Martensite"], -1)
            case 5:
                file1 = get_file_path("Select First XDMF")
                file2 = get_file_path("Select Second XDMF")
                if file1 and file2:
                    ResultPlotting([file1, file2], "Composition/C")

            case 6:
                interactive_datastream_plotter()
            case 7:
                plotting_over_axis()
            case 8:
                checkDatabase()
            case 9:
                MSS.optimization()
            case 10:
                Test()
            case 0:
                print("Exiting...")
                break

            case _:
                print("Option not recognized. Try again.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Goodbye!")
        sys.exit(0)
    print("End")
    sys.exit(0)
