from HelpFile import read_input, checkinput, adjustinputcache
import numpy as np
import meshio

def createMesh():
    if checkinput('Mesh'):
        print('Using precalculated mesh')
        meshdata = meshio.read("Cachefiles/Datastream.xdmf")
        meshio.write("Resultfiles/Datastream.xdmf",
                     meshio.Mesh(points=meshdata.points,
                                 cells={"triangle": meshdata.get_cells_type("triangle")}))
        return
    print('Mesh module')
    data = read_input()
    if data['Programs']["Meshing"] == "Gmsh":
        gmshmodule()
    elif data['Programs']["Meshing"] == "PyGmsh":
        pygmshmodule()
    elif data['Programs']["Meshing"] == "MeshPy":
        meshpymodule()
    else:
        raise KeyError('Meshprogram not implemented in mesh module')
    adjustinputcache('Mesh')
def pygmshmodule():
    print('Meshing with PyGmsh')
    import pygmsh
    import meshio
    import gmsh
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


    meshio.write("Resultfiles/Mesh.nas", mesh)
    #meshio.write("Resultfiles/Mesh.vtk", mesh)
    #meshio.write("Resultfiles/Mesh.msh", mesh)

    # Adding mesh data to datastream
    meshdata = meshio.read("Resultfiles/Mesh.nas")
    meshio.write("Resultfiles/Datastream.xdmf",
                 meshio.Mesh(points=meshdata.points,
                             cells={"triangle": meshdata.get_cells_type("triangle")}))
def gmshmodule():
    print('Meshing with Gmsh')
    import gmsh
    data = read_input()
    r = data['Geometry']['radius']
    tmpgeo = [data['Geometry']['meshscaling'] ** i for i in range(data['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(data['Geometry']['nodes'] - 1)])
    rnodes = r * tmpgeo / np.max(tmpgeo)
    lc = rnodes[-1]-rnodes[-2]
    print(lc)
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("QuarterCirc")
    gdim = 2
    gmsh.option.setNumber("Geometry.Tolerance", 1.E-6)
    gmsh.model.occ.addPoint(0, 0, 0, 1)
    gmsh.model.occ.addPoint(r*np.cos(np.pi/12), r*np.sin(np.pi/12), 0, lc, 2)
    gmsh.model.occ.addPoint(r, 0, 0, lc, 3)
    #gmsh.model.occ.addPoint(r, 0, 0, 3)
    #gmsh.model.occ.addLine(1, 2, 1)
    #gmsh.model.occ.addLine(3, 1, 2)
    #gmsh.model.occ.addCircleArc(2, 1, 3, 3)
    #gmsh.model.occ.addCurveLoop([1, 3, 2], 4)

    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.addLine(3, 1, 2)
    gmsh.model.occ.addCircleArc(2, 1, 3, 3)
    gmsh.model.occ.addCurveLoop([1, 3, 2], 4)

    gmsh.model.occ.addPlaneSurface([4], 1)
    gmsh.model.occ.synchronize()
    #gmsh.model.addPhysicalGroup(1, [1], 1, 'Side')
    #gmsh.model.addPhysicalGroup(1, [2], 2, 'Bottom')
    #gmsh.model.addPhysicalGroup(1, [3], 3, 'Circumference')
    gmsh.model.addPhysicalGroup(gdim, [1], 4, 'Sphere')
    #gmsh.model.mesh.set_size([gdim], 1.E-4)
    #gmsh.model.mesh.set_size_from_boundary(2, 1, 2)
    gmsh.model.mesh.set_transfinite_curve(1, data['Geometry']['nodes'], 'Progression', data['Geometry']['meshscaling'])
    gmsh.model.mesh.set_transfinite_curve(2, data['Geometry']['nodes'], 'Progression', -data['Geometry']['meshscaling'])
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 2*data['Geometry']['radius']/100)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3*data['Geometry']['radius']/100)
    #gmsh.model.mesh.set_order(2)
    #gmsh.option.setNumber('Mesh.ElementOrder', 2)
    #gmsh.model.mesh.optimize()
    gmsh.model.mesh.generate(gdim)
    # ----------------------
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.vtk")
    gmsh.finalize()

    meshdata = meshio.read("Resultfiles/Mesh.msh")
    # Creating datastream from mesh
    meshio.write("Resultfiles/Datastream.xdmf",
                 meshio.Mesh(points=meshdata.points,
                             cells={"triangle": meshdata.get_cells_type("triangle")}))
    data = meshio.read("Resultfiles/Datastream.xdmf")
    # Writing nas file  with meshio
    mesh = meshio.read("Resultfiles/Mesh.msh")
    meshio.write("Resultfiles/Mesh.nas", mesh)

def meshpymodule():
    import meshpy
    import meshpy.triangle as triangle
    import meshpy.tet as tet
    def round_trip_connect(start, end):
        return [(i, i + 1) for i in range(start, end)] + [(end, start)]
    data = read_input()
    r = data['Geometry']['radius']
    lc = data['Geometry']['radius'] * data['Geometry']['meshscaling'] ** (data['Geometry']['nodes'] - 1) / np.sum(
        data['Geometry']['meshscaling'] **
        np.linspace(0, data['Geometry']['nodes'] - 1, data['Geometry']['nodes']))
    tmpgeo = [data['Geometry']['meshscaling']**i for i in range(data['Geometry']['nodes'] - 1)]
    tmpgeo = np.array([np.sum(tmpgeo[0:i]) for i in range(data['Geometry']['nodes'] - 1)])
    rnodes = r*tmpgeo/np.max(tmpgeo)

    #p1 = (0., 0.)
    #p2 = (r * np.cos(np.pi / 6), r * np.sin(np.pi / 6))
    #p3 = (r, 0.)
    p2 = [(rnodes[::-1][i]*np.cos(np.pi / 6), rnodes[::-1][i]*np.sin(np.pi / 6)) for i in range(len(rnodes)-1)]
    p3 = [(rnodes[i], 0.) for i in range(len(rnodes)-1)]
    print(p2[0])
    print(p3[0])
    numel = 30/360*2*r*np.pi/lc
    print(int(numel))
    points = p2+p3

    #
    circ_start = len(points)
    points.extend(
        (r * np.cos(angle), r * np.sin(angle))
        for angle in np.linspace(0, np.pi / 6, int(numel), endpoint=False)
    )
    facets = round_trip_connect(0, len(points)-1)
    #info = tet.MeshInfo()
    #info.set_points(points)
    #info.set_facets(facets)
    #mesh = tet.build()

    #facets.extend(round_trip_connect(circ_start, len(points) - 1))
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)

    mesh = triangle.build(info)
    triangle.write_gnuplot_mesh("Resultfiles/Test.dat", mesh)
    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)

    import matplotlib.pyplot as pt

    pt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)
    pt.show()
