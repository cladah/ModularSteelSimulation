import gmsh
from Framework.HelpFile import read_geninput
from Framework.Datastream_file import createdatastream
import numpy as np
import meshio
import os
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile
import pyvista
import dolfinx
def plot_mesh(mesh: dolfinx.mesh.Mesh, values = None):
    """
    Given a DOLFINx mesh, create a `pyvista.UnstructuredGrid`,
    and plot it and the mesh nodes
    """
    plotter = pyvista.Plotter()
    V_linear = dolfinx.fem.functionspace(mesh, ("Lagrange", 1))
    linear_grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(V_linear))
    if mesh.geometry.cmap.degree > 1:
        ugrid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(mesh))
        if values is not None:
            ugrid.cell_data["Marker"] = values
        plotter.add_mesh(ugrid, style="points", color="b", point_size=10)
        ugrid = ugrid.tessellate()
        plotter.add_mesh(ugrid, show_edges=False)
        plotter.add_mesh(linear_grid,style="wireframe", color="black")

    else:
        if values is not None:
            linear_grid.cell_data["Marker"] = values
        plotter.add_mesh(linear_grid,show_edges=True)
    plotter.show_axes()
    plotter.view_xy()
    if not pyvista.OFF_SCREEN:
        plotter.show()
def gmsh1D():
    print('Remeshing 1D')

    data = read_input()
    r = data['Geometry']['radius']
    tmpgeo = [data['Geometry']['meshscaling'] ** i for i in range(data['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(data['Geometry']['nodes'] - 1)])
    rnodes = r * tmpgeo / np.max(tmpgeo)
    lc = rnodes[-1] - rnodes[-2]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()

    gmsh.clear()
    gmsh.model.add("Line")
    gdim = 1
    gmsh.option.setNumber("Geometry.Tolerance", 1.E-6)

    # Adding three points
    gmsh.model.occ.addPoint(0, 0, 0, 1)
    gmsh.model.occ.addPoint(r, 0, 0, lc, 2)

    # Adding lines and a curveloop between points
    gmsh.model.occ.addLine(1, 2, 1)

    # Adding surface between lines
    gmsh.model.occ.synchronize()

    # Physical group for surface
    gmsh.model.addPhysicalGroup(1, [1], 1, 'Radius')

    # Geometric scaling of mesh
    gmsh.model.mesh.set_transfinite_curve(1, data['Geometry']['nodes'], 'Progression',
                                          data['Geometry']['meshscaling'])

    # Defining element order
    # gmsh.model.mesh.set_order(2)
    # gmsh.model.mesh.recombine()

    # Generating mesh
    gmsh.model.mesh.generate(gdim)
    gmsh.model.mesh.setOrder(1)
    # ----------------------
    if not os.path.exists("Resultfiles"):
        os.makedirs("Resultfiles")
    gmsh.write("Resultfiles/FNXMesh.msh")
    print(*gmsh.logger.get(), sep="\n")
    print(gmsh.model.mesh.get_elements())
    gmsh.finalize()

def four_point_bend(parent):
    ginput = read_geninput()
    gdim = ginput["Geometry"]["dim"]
    if gdim == 2:
        fourPointBend_2D(parent)
    elif gdim == 3:
        fourPointBend_3D(parent)
    else:
        raise KeyError("")

def fourPointBend_3D(parent):
    ginput = read_geninput()
    width = ginput["Geometry"]["width"]
    height = ginput["Geometry"]["height"]
    thick = ginput["Geometry"]["thickness"]
    width_n = int(ginput["Geometry"]["nodesWidth"])
    thick_n = int(ginput["Geometry"]["nodesThickness"])
    height_n = int(ginput["Geometry"]["nodesHeight"])

    tmpgeo = [ginput['Geometry']['meshscaling'] ** i for i in range(ginput['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(ginput['Geometry']['nodes'] - 1)])
    heightnodes = height * tmpgeo / np.max(tmpgeo)
    lc = heightnodes[-1] - heightnodes[-2]

    gmsh.initialize()
    gmsh.logger.start()
    gmsh.model.add("Mesh_3D")
    gdim = 3

    # --- 2D Base Geometry (Same as your 2D code) ---
    p1 = gmsh.model.geo.addPoint(0 , 0, 0)
    p2 = gmsh.model.geo.addPoint(width / 3, 0, 0)
    p3 = gmsh.model.geo.addPoint(width / 2, 0, 0)
    p4 = gmsh.model.geo.addPoint(width / 2, height, 0)
    p5 = gmsh.model.geo.addPoint(width / 3, height, 0)
    p6 = gmsh.model.geo.addPoint(0 , height, 0)
    p7 = gmsh.model.geo.addPoint(-0.01, 0, 0)
    p8 = gmsh.model.geo.addPoint(-0.01, height, 0)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p5)
    l3 = gmsh.model.geo.addLine(p5, p6)
    l4 = gmsh.model.geo.addLine(p6, p1)
    l5 = gmsh.model.geo.addLine(p2, p3)
    l6 = gmsh.model.geo.addLine(p3, p4)
    l7 = gmsh.model.geo.addLine(p4, p5)
    l8 = gmsh.model.geo.addLine(p7, p1)
    l9 = gmsh.model.geo.addLine(p6, p8)
    l10 = gmsh.model.geo.addLine(p8, p7)

    s1 = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])])
    s2 = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop([l5, l6, l7, -l2])])
    s3 = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop([l8, -l4, l9, l10])])

    # --- Transfinite Constraints for the 2D Face ---
    for l in [l8, l9]:
        gmsh.model.geo.mesh.setTransfiniteCurve(l, int(4))
    for l in [l1, l3]:
        gmsh.model.geo.mesh.setTransfiniteCurve(l, int(width_n))
    for l in [l5, l7]:
        gmsh.model.geo.mesh.setTransfiniteCurve(l, int(width_n / 2))
    for l in [l2, l4, l6, l10]:
        gmsh.model.geo.mesh.setTransfiniteCurve(l, height_n, 'Bump', 0.1)

    gmsh.model.geo.mesh.setTransfiniteSurface(s1)
    gmsh.model.geo.mesh.setTransfiniteSurface(s2)
    gmsh.model.geo.mesh.setTransfiniteSurface(s3)
    gmsh.model.geo.mesh.setRecombine(2, s1)
    gmsh.model.geo.mesh.setRecombine(2, s2)
    gmsh.model.geo.mesh.setRecombine(2, s3)

    ext1 = gmsh.model.geo.extrude([(2, s1), (2, s2), (2, s3)], 0, 0, thick/2, numElements=[thick_n - 1], recombine=True)
    gmsh.model.geo.synchronize()
    # Extract the volume IDs from the extrusion results
    volumes = [item[1] for item in ext1 if item[0] == 3]
    gmsh.model.addPhysicalGroup(3, volumes, 4, "Grid")

    parent.updateprogress(0.5)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(gdim)
    #gmsh.model.mesh.setOrder(2)
    #gmsh.fltk.run()
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.vtk")
    print(*gmsh.logger.get(), sep="\n")
    parent.updateprogress(0.8)

    mesh = meshio.read("Resultfiles/Mesh.msh")
    #meshio.write("Resultfiles/Mesh.nas", mesh)
    domain, cell_tags, facet_tags = gmshio.model_to_mesh(
        gmsh.model, MPI.COMM_WORLD, 0, gdim=gdim
    )
    domain.name = "Grid"
    with XDMFFile(domain.comm, "Resultfiles/Mesh.xdmf", "w") as file:
        file.write_mesh(domain)

    createdatastream(mesh)
    #gmsh.fltk.run()
    gmsh.finalize()


    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()

    # Plot with PyVista
    #plotter = pv.Plotter()
    #plotter.add_mesh(grid, show_edges=True)
    #plotter.show(
def fourPointBend_2D(parent):
    ginput = read_geninput()
    width = ginput["Geometry"]["width"]
    height = ginput["Geometry"]["height"]
    width_n = int(ginput["Geometry"]["nodesWidth"])
    print(width_n/3)
    height_n = int(ginput["Geometry"]["nodesHeight"])

    tmpgeo = [ginput['Geometry']['meshscaling'] ** i for i in range(ginput['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(ginput['Geometry']['nodes'] - 1)])
    heightnodes = height * tmpgeo / np.max(tmpgeo)
    lc = heightnodes[-1] - heightnodes[-2]

    gmsh.initialize()
    #gmsh.option.setNumber("General.Terminal", 0)

    gmsh.logger.start()

    gmsh.clear()
    gmsh.model.add("Mesh")
    gdim = 2
    #gmsh.option.setNumber("Geometry.Tolerance", 1.E-6)
    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint(width/3, 0, 0)
    p3 = gmsh.model.geo.addPoint(2*width/3, 0, 0)
    p4 = gmsh.model.geo.addPoint(width, 0, 0)
    p5 = gmsh.model.geo.addPoint(width, height, 0)
    p6 = gmsh.model.geo.addPoint(2*width/3, height, 0)
    p7 = gmsh.model.geo.addPoint(width/3, height, 0)
    p8 = gmsh.model.geo.addPoint(0, height, 0)

    l1 = gmsh.model.geo.addLine(p1,p2)
    l2 = gmsh.model.geo.addLine(p2,p7)
    l3 = gmsh.model.geo.addLine(p7,p8)
    l4 = gmsh.model.geo.addLine(p8,p1)
    l5 = gmsh.model.geo.addLine(p2,p3)
    l6 = gmsh.model.geo.addLine(p3,p6)
    l7 = gmsh.model.geo.addLine(p6,p7)
    l8 = gmsh.model.geo.addLine(p3,p4)
    l9 = gmsh.model.geo.addLine(p4,p5)
    l10 = gmsh.model.geo.addLine(p5,p6)
    cl1 = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4])
    cl2 = gmsh.model.geo.addCurveLoop([l5,l6,l7,-l2])
    cl3 = gmsh.model.geo.addCurveLoop([l8,l9,l10,-l6])
    s1 = gmsh.model.geo.addPlaneSurface([cl1])
    s2 = gmsh.model.geo.addPlaneSurface([cl2])
    s3 = gmsh.model.geo.addPlaneSurface([cl3])
    parent.updateprogress(0.5)

    # Geometric scaling of mesh

    gmsh.model.geo.mesh.setTransfiniteCurve(l1, int(width_n/3))
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, height_n, 'Bump', 0.1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, int(width_n/3))
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, height_n, 'Bump', 0.1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l5, int(width_n/3))
    gmsh.model.geo.mesh.setTransfiniteCurve(l6, height_n, 'Bump', 0.1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l7, int(width_n/3))
    gmsh.model.geo.mesh.setTransfiniteCurve(l8, int(width_n/3))
    gmsh.model.geo.mesh.setTransfiniteCurve(l9, height_n, 'Bump', 0.1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l10, int(width_n/3))

    gmsh.model.geo.mesh.setTransfiniteSurface(s1)
    gmsh.model.geo.mesh.setTransfiniteSurface(s2)
    gmsh.model.geo.mesh.setTransfiniteSurface(s3)
    gmsh.model.geo.mesh.setRecombine(2, s1)  # Quadrilateral elements
    gmsh.model.geo.mesh.setRecombine(2, s2)
    gmsh.model.geo.mesh.setRecombine(2, s3)
    #compound_surface = gmsh.model.geo.addSurfaceLoop([s1, s2, s3])
    #s_combined = gmsh.model.geo.addPlaneSurface([compound_surface])
    #gmsh.model.geo.mesh.setTransfiniteSurface(s_combined)
    #gmsh.model.geo.mesh.setRecombine(2, s_combined)
    #gmsh.model.addPhysicalGroup(2, [s_combined], 4, "Grid")
    gmsh.model.addPhysicalGroup(2, [s1,s2,s3], 4, "Grid")

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(gdim)
    #gmsh.model.mesh.setOrder(2)
    #gmsh.fltk.run()
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.vtk")
    print(*gmsh.logger.get(), sep="\n")
    parent.updateprogress(0.8)

    mesh = meshio.read("Resultfiles/Mesh.msh")
    meshio.write("Resultfiles/Mesh.nas", mesh)
    domain, cell_tags, facet_tags = gmshio.model_to_mesh(
        gmsh.model, MPI.COMM_WORLD, 0, gdim=2
    )
    domain.name = "Grid"
    with XDMFFile(domain.comm, "Resultfiles/Mesh.xdmf", "w") as file:
        file.write_mesh(domain)

    createdatastream(mesh)
    gmsh.finalize()


    with meshio.xdmf.TimeSeriesReader("Datastream.xdmf") as reader:
        points, cells = reader.read_points_cells()

    # Plot with PyVista
    #plotter = pv.Plotter()
    #plotter.add_mesh(grid, show_edges=True)
    #plotter.show()

def cylinder(parent):
    ginput = read_geninput()
    r = ginput['Geometry']['radius']
    tmpgeo = [ginput['Geometry']['meshscaling'] ** i for i in range(ginput['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(ginput['Geometry']['nodes'] - 1)])
    rnodes = r * tmpgeo / np.max(tmpgeo)
    lc = rnodes[-1] - rnodes[-2]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()

    gmsh.clear()
    gmsh.model.add("Grid")
    gdim = 2
    gmsh.option.setNumber("Geometry.Tolerance", 1.E-6)

    # Adding three points
    gmsh.model.occ.addPoint(0, 0, 0, 1)
    gmsh.model.occ.addPoint(r * np.cos(np.pi / 24), r * np.sin(np.pi / 24), 0, lc, 2)
    gmsh.model.occ.addPoint(r, 0, 0, lc, 3)

    # Adding lines and a curveloop between points
    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.addLine(3, 1, 2)
    gmsh.model.occ.addCircleArc(2, 1, 3, 3)
    gmsh.model.occ.addCurveLoop([1, 3, 2], 4)

    # Adding surface between lines
    gmsh.model.occ.addPlaneSurface([4], 1)
    gmsh.model.occ.synchronize()
    parent.updateprogress(0.5)

    # Physical group for surface
    gmsh.model.addPhysicalGroup(gdim, [1], 4, 'Sphere')

    # Geometric scaling of mesh
    gmsh.model.mesh.set_transfinite_curve(1, ginput['Geometry']['nodes'], 'Progression',
                                          ginput['Geometry']['meshscaling'])
    gmsh.model.mesh.set_transfinite_curve(2, ginput['Geometry']['nodes'], 'Progression',
                                          -ginput['Geometry']['meshscaling'])


    gmsh.model.mesh.generate(gdim)

    gmsh.model.mesh.setOrder(2)

    # ----------------------
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.vtk")
    # print(gmsh.model.mesh.setOrder())
    print(*gmsh.logger.get(), sep="\n")

    gmsh.finalize()
    parent.updateprogress(0.8)
    mesh = meshio.read("Resultfiles/Mesh.msh")
    meshio.write("Resultfiles/Mesh.nas", mesh)

    len(mesh.get_cells_type("triangle6"))
    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(points=mesh.points[:, :2], cells={"triangle6": mesh.get_cells_type("triangle6")})
        writer.write_data(0, cell_data={})
    # Writing nas file with meshio (Needed for Comsol FEM solver)

    mesh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, gdim=2)
    mesh.name = "Grid"
    createdatastream(mesh)

    parent.updateprogress(1.0)



def gmshsolver(parent):
    print('Meshing with Gmsh')
    ginput = read_geninput()
    geometries = ["2Daxisym", "4PointBend"]
    if ginput["Geometry"]["Type"] == "2Daxisym":
        cylinder(parent)
    elif ginput["Geometry"]["Type"] == "4PointBend":
        four_point_bend(parent)
    else:
        raise KeyError("Geometry not implemented. ")


def pygmshsolver(parent):
    print('Meshing with PyGmsh')
    data = read_input()
    r = data['Geometry']['radius']
    tmpgeo = [data['Geometry']['meshscaling'] ** i for i in range(data['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(data['Geometry']['nodes'] - 1)])
    rnodes = r * tmpgeo / np.max(tmpgeo)
    lc = rnodes[-1] - rnodes[-2]

    with pygmsh.geo.Geometry() as geom:
        p0 = geom.add_point((0.0, 0.0, 0.0))
        p1 = geom.add_point((r, 0.0, 0.0), mesh_size=lc)
        p2 = geom.add_point((r*np.cos(np.pi/6), r*np.sin(np.pi/6), 0.0),mesh_size =lc)
        l1 = geom.add_line(p0, p1)
        l2 = geom.add_line(p2, p0)
        l3 = geom.add_circle_arc(p1, p0, p2)
        loop1 = geom.add_curve_loop((l1, l3, l2))
        rs0 = geom.add_surface(loop1)
        geom.synchronize()
        # Adding Geometric scaling of the mesh
        geom.set_transfinite_curve(l1,data["Geometry"]["nodes"],"Progression",data["Geometry"]["meshscaling"])
        geom.set_transfinite_curve(l2, data["Geometry"]["nodes"], "Progression", -data["Geometry"]["meshscaling"])
        geom.add_physical([rs0], "Volume")

        mesh = geom.generate_mesh(dim=2)
    mesh = mesh.prune_z_0()
    print(mesh)
    meshio.write("Resultfiles/Mesh.nas", mesh)
    #meshio.write("Resultfiles/Mesh.vtk", mesh)
    #meshio.write("Resultfiles/Mesh.msh", mesh)

    # Adding mesh data to datastream
    meshdata = meshio.read("Resultfiles/Mesh.nas")
    meshio.write("Resultfiles/Datastream.xdmf",
                 meshio.Mesh(points=meshdata.points,
                             cells={"triangle": meshdata.get_cells_type("triangle")}))