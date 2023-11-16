import mph
import os

def runComsol():
    directory = os.getcwd()
    meshdirec = directory + '\\Resultfiles\\Mesh.nas'
    savedirec = directory + '\\Resultfiles'
    client = mph.start()
    pymodel = client.load("Resultfiles/Comsolmodel.mph")
    model = pymodel.java

    model.func("int1").setIndex("table", 100., 0, 0)
    model.func("int1").setIndex("table", 100., 1, 0)

    model.save('Resultfiles/Comsolmodel')
def setupComsol():
    directory = os.getcwd()
    meshdirec = directory + '\\Resultfiles\\Mesh.nas'
    savedirec = directory + '\\Resultfiles'
    client = mph.start()
    pymodel = client.create()
    model = pymodel.java
    model.modelNode().create("comp1", True)

    # --------------- Setting up geometry and mesh ------------------#
    model.component("comp1").geom().create("geom1", 2)
    model.component("comp1").geom("geom1").axisymmetric(True)
    model.component("comp1").mesh().create("mesh1")
    model.component("comp1").geom("geom1").run()
    model.component("comp1").mesh("mesh1").create("imp1", "Import")
    model.component("comp1").mesh("mesh1").feature("imp1").set("source", "nastran")
    model.component("comp1").mesh("mesh1").feature("imp1").set("filename", meshdirec)
    model.component("comp1").mesh("mesh1").feature("imp1").importData()
    model.component("comp1").mesh("mesh1").run()

    # --------------- Setting up parameters ------------------#
    model.func().create("EM", "Interpolation");
    model.func().create("EA", "Interpolation");
    model.func().create("EP", "Interpolation");
    model.func().create("EB", "Interpolation");

    # --------------- Setting up material ------------------#
    model.component("comp1").material().create("159", "Common")
    model.component("comp1").material("159").propertyGroup("def").set("youngsmodulus", "210E9")
    model.component("comp1").material("159").propertyGroup("def").set("poissonsratio", "0.3")
    model.component("comp1").material("159").propertyGroup("def").set("thermalconductivity", "45")
    model.component("comp1").material("159").propertyGroup("def").set("density", "7800")
    model.component("comp1").material("159").propertyGroup("def").set("heatcapacity", "470")
    model.component("comp1").material("159").selection().set(1)

    # --------------- Setting up phase transformation ------------------#
    model.func().create("an1", "Analytic")
    model.func("an1").set("funcname", "Ms")
    model.func("an1").set("expr", "0.011+0.01*carb[1/K]")
    model.func("an1").set("args", "carb")
    model.func("an1").set("expr", "200+100*carb[degC]")

    # --------------- Adding physics ------------------#
    model.component("comp1").physics().create("solid", "SolidMechanics", "geom1")
    model.component("comp1").physics().create("ht", "HeatTransfer", "geom1")
    model.component("comp1").physics().create("audc", "AusteniteDecomposition", "geom1")

    # --------------- Adjusting properties
    model.component("comp1").physics("ht").prop("ShapeProperty").set("order_temperature", "1")
    #
    model.component("comp1").physics("audc").prop("SolidMechanics").set("trip", "1")
    model.component("comp1").physics("audc").prop("SolidMechanics").set("plasticity", "0")
    model.component("comp1").physics("audc").prop("SolidMechanics").set("dilstrain", "1")
    model.component("comp1").physics("audc").prop("HeatTransfer").set("latentheat", "1")
    model.component("comp1").physics("audc").prop("ShapeProperty").set("order_straindiscr_disc", "2")

    # --------------- Adding multiphysics ------------------#
    model.component("comp1").multiphysics().create("lht1", "PhaseTransformationLatentHeat", 2)
    model.component("comp1").multiphysics("lht1").set("Metphase_physics", "audc")
    model.component("comp1").multiphysics("lht1").set("HeatTransfer_physics", "ht")
    model.component("comp1").multiphysics("lht1").selection().all()
    model.component("comp1").multiphysics().create("ptstr1", "PhaseTransformationStrain", 2)
    model.component("comp1").multiphysics("ptstr1").set("Metphase_physics", "audc")
    model.component("comp1").multiphysics("ptstr1").set("SolidMechanics_physics", "solid")
    model.component("comp1").multiphysics("ptstr1").selection().all()

    # --------------- Adding boundary conditions ------------------#
    model.component("comp1").physics("solid").create("symp1", "SymmetryPlane", 1)
    model.component("comp1").physics("solid").feature("symp1").selection().set(2)

    model.component("comp1").physics("ht").create("symp1", "Symmetry", 1)
    model.component("comp1").physics("ht").feature("symp1").selection().set(2)
    model.component("comp1").physics("ht").create("hf1", "HeatFluxBoundary", 1)
    model.component("comp1").physics("ht").feature("hf1").selection().set(3)
    model.component("comp1").physics("ht").feature("hf1").set("q0_input", 1234.)

    #model.component("comp1").physics("ht").feature("hf1").set("q0_input", 12345)  # Add the correct temperature inflow

    #model.component("comp1").physics("ht").create("sym1", "Symmetry", 1)
    #model.component("comp1").physics("ht").feature("sym1").selection().set(2)

    # --------------- Creating study ------------------#
    model.component("comp1").physics("audc").active(False);
    model.component("comp1").multiphysics("lht1").active(False);
    model.component("comp1").multiphysics("ptstr1").active(False);

    model.save('Resultfiles/Comsolmodel')
    return
