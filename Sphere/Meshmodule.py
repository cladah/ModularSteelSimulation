import meshio

from HelpFile import read_input, checkinput, adjustinputcache
import numpy as np
import meshio

def createMesh():
    if checkinput('Mesh'):
        print('Using precalculated mesh')
        return
    print('Mesh module')
    data = read_input()
    if data['Programs']["Meshing"] == "Gmsh":
        gmshmodule()
        #pygmshmodule()
    else:
        raise KeyError('Meshprogram not implemented in mesh module')
    adjustinputcache('Mesh')
def pygmshmodule():
    print('Meshing with PyGmsh')
    import pygmsh
    import meshio
    data = read_input()
    r = data['Geometry']['radius']
    lc = data['Geometry']['radius'] * data['Geometry']['meshscaling'] ** (data['Geometry']['nodes'] - 1) / np.sum(
        data['Geometry']['meshscaling'] **
        np.linspace(0, data['Geometry']['nodes'] - 1, data['Geometry']['nodes']))

    with pygmsh.geo.Geometry() as geom:
        p0 = geom.add_point((0.0, 0.0, 0.0))
        p1 = geom.add_point((r, 0.0, 0.0))
        p2 = geom.add_point((r*np.cos(np.pi/6), r*np.sin(np.pi/6), 0.0),lc)
        l1 = geom.add_line(p0, p1)
        l2 = geom.add_line(p2, p0)
        l3 = geom.add_circle_arc(p1, p0, p2)
        loop1 = geom.add_curve_loop((l1, l3, l2))
        rs0 = geom.add_surface(loop1)
        geom.set_transfinite_curve(l1,data["Geometry"]["nodes"],"Progression",data["Geometry"]["meshscaling"])
        geom.set_transfinite_curve(l2, data["Geometry"]["nodes"], "Progression", -data["Geometry"]["meshscaling"])

        #geom.set_transfinite_curve([l1], progression=data['Geometry']['meshscaling'])
        #geom.set_transfinite_curve([l2])
        #geom.set_transfinite_surface(rs0)
        #geom.add_circle_arc([r, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, r, 0.0])
        #geom.add_circle_arc([r, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0,r,0.0])
        #geom.add_circle([0.0, 0.0], 1.0, mesh_size=0.2)
        mesh = geom.generate_mesh()
    #meshio.write("test.nas",mesh)
def gmshmodule():
    print('Meshing with Gmsh')
    import gmsh
    data = read_input()
    lc = data['Geometry']['radius'] * data['Geometry']['meshscaling'] ** (data['Geometry']['nodes'] - 1) / np.sum(data['Geometry']['meshscaling'] **
        np.linspace(0, data['Geometry']['nodes'] - 1, data['Geometry']['nodes']))

    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("QuarterCirc")
    gdim = 2
    r = data['Geometry']['radius']
    gmsh.model.occ.addPoint(0, 0, 0, 1)
    gmsh.model.occ.addPoint(r*np.cos(np.pi/6), r*np.sin(np.pi/6), 0, lc, 2)
    gmsh.model.occ.addPoint(r, 0, 0, lc, 3)
    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.addLine(1, 3, 2)
    gmsh.model.occ.addCircleArc(2, 1, 3, 3)
    gmsh.model.occ.addCurveLoop([1, 2, 3], 4)
    gmsh.model.occ.addPlaneSurface([4], 1)
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(1, [1], 1, 'Bottom')
    gmsh.model.addPhysicalGroup(1, [2], 2, 'Side')
    gmsh.model.addPhysicalGroup(1, [3], 3, 'Circumference')
    gmsh.model.addPhysicalGroup(gdim, [1], 4, 'Sphere')
    gmsh.model.mesh.set_transfinite_curve(1, data['Geometry']['nodes'], 'Progression', data['Geometry']['meshscaling'])
    gmsh.model.mesh.set_transfinite_curve(2, data['Geometry']['nodes'], 'Progression', data['Geometry']['meshscaling'])
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 2*data['Geometry']['radius']/100)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3*data['Geometry']['radius']/100)
    #gmsh.model.mesh.set_order(2)
    gmsh.option.setNumber('Mesh.ElementOrder', 2)
    gmsh.model.mesh.optimize()
    gmsh.model.mesh.generate(gdim)
    # ----------------------
    gmsh.write("Resultfiles/Mesh.msh")
    gmsh.write("Resultfiles/Mesh.nas")
    gmsh.write("Resultfiles/Mesh.vtk")
    gmsh.finalize()

    meshdata = meshio.read("Resultfiles/Mesh.msh")
    meshio.write("Resultfiles/Datastream.xdmf",
                 meshio.Mesh(points=meshdata.points,
                             cells={"triangle": meshdata.get_cells_type("triangle")}))
                             #point_data={"T": 0.*np.arange(len(meshdata.points)),"C": 0.*np.arange(len(meshdata.points))},
                             #cell_data={"phiM": [0.*np.arange(len(meshdata.cells[3]))]}))
    data = meshio.read("Resultfiles/Datastream.xdmf")
