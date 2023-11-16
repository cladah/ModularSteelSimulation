import mph
import os


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
    phases = ["Austenite", "Ferrite", "Perlite", "Bainite", "Martensite"]
    for phase in phases:
        ph = phase[0]
        model.nodeGroup().create(phase, "GlobalDefinitions")
        model.func().create("E_"+ph, "Interpolation")
        model.nodeGroup(phase).add("func", "E_"+ph)
        model.func().create("Sy_"+ph, "Interpolation")
        model.nodeGroup(phase).add("func", "Sy_"+ph)
        model.func().create("Cp_"+ph, "Interpolation")
        model.nodeGroup(phase).add("func", "Cp_"+ph)
        model.func().create("k_"+ph, "Interpolation")
        model.nodeGroup(phase).add("func", "k_"+ph)
        model.func().create("alpha_"+ph, "Interpolation")
        model.nodeGroup(phase).add("func", "alpha_"+ph)
        if phase == "Austenite":
            pass
        elif phase == "Martensite":
            model.func().create("Ms", "Interpolation")
            model.nodeGroup("Martensite").add("func", "Ms")
            model.func().create("beta", "Interpolation")
            model.nodeGroup("Martensite").add("func", "beta")
        else:
            model.func().create("Tau_"+ph, "Interpolation")
            model.nodeGroup(phase).add("func", "Tau_"+ph)
            model.func().create("n_"+ph, "Interpolation")
            model.nodeGroup(phase).add("func", "n_"+ph)

    # --------------- Setting up material ------------------#
    for phase in phases:
        ph = phase[0]
        model.material().create(phase, "Common", "")
        model.material(phase).propertyGroup("def").set("heatcapacity", "Cp_" + ph)
        model.material(phase).propertyGroup("def").set("thermalconductivity", "k_" + ph)
        model.material(phase).propertyGroup().create("ThermalExpansion", "Thermal_expansion")
        model.material(phase).propertyGroup("ThermalExpansion").set("thermalexpansioncoefficient", "alpha_" + ph)
        model.material(phase).propertyGroup().create("Enu", "Young's_modulus_and_Poisson's_ratio")
        model.material(phase).propertyGroup("Enu").set("E", "E_" + ph)
        model.material(phase).propertyGroup("Enu").set("density", "7800")
        model.material(phase).propertyGroup("Enu").set("nu", "0.3")
        model.material(phase).propertyGroup().create("ElastoplasticModel", "Elastoplastic material model")
        model.material(phase).propertyGroup("ElastoplasticModel").set("sigmags", "Sy_" + ph)

    # --------------- Adding physics ------------------#
    model.component("comp1").physics().create("solid", "SolidMechanics", "geom1")
    model.component("comp1").physics().create("ht", "HeatTransfer", "geom1")
    model.component("comp1").physics().create("audc", "AusteniteDecomposition", "geom1")

    model.component("comp1").physics("audc").prop("MaterialProperties").runCommand("makecompoundmaterial")
    model.component("comp1").material("audcmat").propertyGroup("def").set("thermalexpansioncoefficient", "comp1.audc.alpha_iso")


    # Austenite
    model.component("comp1").physics("audc").feature("phase1").set("phaseMaterial", "Austenite")

    # Ferrite
    model.component("comp1").physics("audc").feature("phase2").set("phaseMaterial", "Ferrite")
    model.component("comp1").physics("audc").feature("ptran1").set("ptModel", "JMAK")
    model.component("comp1").physics("audc").feature("ptran1").set("taujmak", "Tau_F")
    model.component("comp1").physics("audc").feature("ptran1").set("njmak", "n_F")
    #model.component("comp1").physics("audc").feature("ptran1").set("trip", True)

    # Perlite
    model.component("comp1").physics("audc").feature("phase3").set("phaseMaterial", "Perlite")
    model.component("comp1").physics("audc").feature("ptran2").set("ptModel", "JMAK")
    model.component("comp1").physics("audc").feature("ptran2").set("taujmak", "Tau_P")
    model.component("comp1").physics("audc").feature("ptran2").set("njmak", "n_P")
    #model.component("comp1").physics("audc").feature("ptran2").set("trip", True)

    # Bainite
    model.component("comp1").physics("audc").feature("phase2").set("phaseMaterial", "Ferrite")
    model.component("comp1").physics("audc").feature("ptran3").set("ptModel", "JMAK")
    model.component("comp1").physics("audc").feature("ptran3").set("taujmak", "Tau_B")
    model.component("comp1").physics("audc").feature("ptran3").set("njmak", "n_B")
    #model.component("comp1").physics("audc").feature("ptran3").set("trip", True)

    # Martensite
    model.component("comp1").physics("audc").feature("ptran4").set("Ms", "Ms")
    model.component("comp1").physics("audc").feature("ptran4").set("beta", "beta")

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
    model.component("comp1").multiphysics().create("te1", "ThermalExpansion", 2)
    model.component("comp1").multiphysics("te1").selection().all()

    # --------------- Adding boundary conditions ------------------#
    model.component("comp1").physics("solid").create("symp1", "SymmetryPlane", 1)
    model.component("comp1").physics("solid").feature("symp1").selection().set(2)

    model.component("comp1").physics("ht").create("symp1", "Symmetry", 1)
    model.component("comp1").physics("ht").feature("symp1").selection().set(2)
    model.component("comp1").physics("ht").create("hf1", "HeatFluxBoundary", 1)
    model.component("comp1").physics("ht").feature("hf1").selection().set(3)
    model.component("comp1").physics("ht").feature("hf1").set("q0_input", -1234.)

    #model.component("comp1").physics("ht").feature("hf1").set("q0_input", 12345)  # Add the correct temperature inflow

    #model.component("comp1").physics("ht").create("sym1", "Symmetry", 1)
    #model.component("comp1").physics("ht").feature("sym1").selection().set(2)




    # --------------- Setting up phase transformation ------------------#
    # model.func().create("an1", "Analytic")
    # model.func("an1").set("funcname", "Ms")
    # model.func("an1").set("expr", "0.011+0.01*carb[1/K]")
    # model.func("an1").set("args", "carb")
    # model.func("an1").set("expr", "200+100*carb[degC]")

    # --------------- Creating study ------------------#
    #model.component("comp1").physics("audc").active(False)
    #model.component("comp1").multiphysics("lht1").active(False)
    #model.component("comp1").multiphysics("ptstr1").active(False)

    model.save('Resultfiles/Comsolmodel')
    return
def setupComsolSolver():
    pass

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

    setupComsolSolver()