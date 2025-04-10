import gmsh
from Framework.HelpFile import read_geninput
import numpy as np
import meshio
import pygmsh
import os

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
def gmshsolver(parent):
    print('Meshing with Gmsh')

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
    gmsh.model.add("QuarterCirc")
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

    # PhysicalGroup if needed
    # gmsh.model.addPhysicalGroup(1, [1], 1, 'Side')
    # gmsh.model.addPhysicalGroup(1, [2], 2, 'Bottom')
    # gmsh.model.addPhysicalGroup(1, [3], 3, 'Circumference')

    # Physical group for surface
    gmsh.model.addPhysicalGroup(gdim, [1], 4, 'Sphere')

    # Geometric scaling of mesh
    gmsh.model.mesh.set_transfinite_curve(1, ginput['Geometry']['nodes'], 'Progression',
                                          ginput['Geometry']['meshscaling'])
    gmsh.model.mesh.set_transfinite_curve(2, ginput['Geometry']['nodes'], 'Progression',
                                          -ginput['Geometry']['meshscaling'])

    # Defining element order

    # gmsh.model.mesh.set_order(2)
    # gmsh.model.mesh.recombine()

    # Generating mesh
    gmsh.model.mesh.generate(gdim)

    #
    gmsh.model.mesh.setOrder(2)
    # ----------------------
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.vtk")
    #print(gmsh.model.mesh.setOrder())
    print(*gmsh.logger.get(), sep="\n")
    gmsh.finalize()
    parent.updateprogress(0.8)
    mesh = meshio.read("Resultfiles/Mesh.msh")
    meshio.write("Resultfiles/Mesh.nas", mesh)


    # -----------
    # reader = vtk.vtkUnstructuredGridReader()
    # reader.SetFileName("Resultfiles/Mesh.vtk")
    # reader.Update()
    #
    # # Get the unstructured grid from the reader
    # mesh = reader.GetOutput()
    #
    # writer = vtk.vtkXdmfWriter()
    # writer.SetFileName("Datastream.xdmf")
    # writer.SetInputData(mesh)
    #
    # # Write XDMF file
    # writer.Write()
    # ---------------------------------


    # Creating datastream from mesh
    # meshio.write("Datastream.xdmf",
    #              meshio.Mesh(points=mesh.points,
    #                          cells={"triangle": mesh.get_cells_type("triangle")}))
    len(mesh.get_cells_type("triangle6"))
    with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
        writer.write_points_cells(points=mesh.points[:, :2], cells={"triangle6": mesh.get_cells_type("triangle6")})
        writer.write_data(0, cell_data={})
    # Writing nas file with meshio (Needed for Comsol FEM solver)

    parent.updateprogress(1.0)

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