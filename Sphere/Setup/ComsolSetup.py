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

    #model.sol().create("sol1");

    #model.sol("sol1").study("std1")
    model.study().create("std1")
    model.study("std1").create("time", "Transient")
    model.study("std1").feature("time").setSolveFor("/physics/solid", True)
    model.study("std1").feature("time").setSolveFor("/physics/ht", True)
    model.study("std1").feature("time").setSolveFor("/physics/audc", False)
    model.study("std1").feature("time").setSolveFor("/multiphysics/lht1", False)
    model.study("std1").feature("time").setSolveFor("/multiphysics/ptstr1", False)
    model.study("std1").setGenPlots(False);
    model.study("std1").setGenConv(False)

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    # model.study("std1").feature("time").set("notlistsolnum", 1);
    model.study("std1").feature("time").set("notsolnum", "auto");
    # model.study("std1").feature("time").set("listsolnum", 1);
    model.study("std1").feature("time").set("solnum", "auto");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").feature("st1").set("study", "std1");
    model.sol("sol1").feature("st1").set("studystep", "time");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").feature("v1").feature("comp1_u").set("scalemethod", "manual");
    model.sol("sol1").feature("v1").feature("comp1_u").set("scaleval", "1e-2*0.14142135623730953");
    model.sol("sol1").feature("v1").set("control", "time");
    model.sol("sol1").create("t1", "Time");
    model.sol("sol1").feature("t1").set("tlist", "range(0,0.1,1)");
    model.sol("sol1").feature("t1").set("plot", "off");
    model.sol("sol1").feature("t1").set("plotgroup", "Default");
    model.sol("sol1").feature("t1").set("plotfreq", "tout");
    model.sol("sol1").feature("t1").set("probesel", "all");
    model.sol("sol1").feature("t1").set("probes", []);
    model.sol("sol1").feature("t1").set("probefreq", "tsteps");
    model.sol("sol1").feature("t1").set("rtol", 0.001);
    model.sol("sol1").feature("t1").set("atolglobalvaluemethod", "factor");
    model.sol("sol1").feature("t1").set("atolmethod", ["comp1_T", "global", "comp1_u", "global"]);
    model.sol("sol1").feature("t1").set("atol", ["comp1_T", "1e-3", "comp1_u", "1e-3"]);
    model.sol("sol1").feature("t1").set("atolvaluemethod", ["comp1_T", "factor", "comp1_u", "factor"]);
    model.sol("sol1").feature("t1").set("atolvaluemethod", ["comp1_T", "factor", "comp1_u", "factor"]);
    model.sol("sol1").feature("t1").set("reacf", True);
    model.sol("sol1").feature("t1").set("storeudot", True);
    model.sol("sol1").feature("t1").set("endtimeinterpolation", True);
    model.sol("sol1").feature("t1").set("timemethod", "genalpha");
    model.sol("sol1").feature("t1").set("estrat", "exclude");
    model.sol("sol1").feature("t1").set("maxorder", 2);
    model.sol("sol1").feature("t1").set("control", "time");
    model.sol("sol1").feature("t1").feature("aDef").set("cachepattern", True);
    model.sol("sol1").feature("t1").create("seDef", "Segregated");
    model.sol("sol1").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("t1").feature("fc1").set("jtech", "once");
    model.sol("sol1").feature("t1").feature("fc1").set("damp", 0.9);
    model.sol("sol1").feature("t1").feature("fc1").set("stabacc", "aacc");
    model.sol("sol1").feature("t1").feature("fc1").set("aaccdim", 5);
    model.sol("sol1").feature("t1").feature("fc1").set("aaccmix", 0.9);
    model.sol("sol1").feature("t1").feature("fc1").set("aaccdelay", 0);
    model.sol("sol1").feature("t1").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("d1").set("linsolver", "mumps");
    model.sol("sol1").feature("t1").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature("d1").label("Direct, heat transfer variables (ht) (merged)");
    model.sol("sol1").feature("t1").create("i1", "Iterative");
    model.sol("sol1").feature("t1").feature("i1").set("linsolver", "gmres");
    model.sol("sol1").feature("t1").feature("i1").set("prefuntype", "left");
    model.sol("sol1").feature("t1").feature("i1").set("itrestart", 50);
    model.sol("sol1").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i1").set("maxlinit", 10000);
    model.sol("sol1").feature("t1").feature("i1").set("nlinnormuse", "on");
    model.sol("sol1").feature("t1").feature("i1").label("AMG, heat transfer variables (ht)");
    model.sol("sol1").feature("t1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("prefun", "saamg");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("mgcycle", "v");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("maxcoarsedof", 50000);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("strconn", 0.01);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("nullspace", "constant");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("usesmooth", False);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("saamgcompwise", True);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").set("loweramg", True);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("so1").set("iter", 2);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("pr").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").create("so1", "SOR");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("so1").set("iter", 2);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("po").feature("so1").set("relax", 0.9);
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").set("linsolver",
                                                                                                 "pardiso");

    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").set("pivotperturb",
                                                                                                 1.0E-13);
    model.sol("sol1").feature("t1").feature("fc1").set("linsolver", "d1");
    model.sol("sol1").feature("t1").feature("fc1").set("jtech", "once");
    model.sol("sol1").feature("t1").feature("fc1").set("damp", 0.9);
    model.sol("sol1").feature("t1").feature("fc1").set("stabacc", "aacc");
    model.sol("sol1").feature("t1").feature("fc1").set("aaccdim", 5);
    model.sol("sol1").feature("t1").feature("fc1").set("aaccmix", 0.9);
    model.sol("sol1").feature("t1").feature("fc1").set("aaccdelay", 0);
    model.sol("sol1").feature("t1").feature().remove("fcDef");
    model.sol("sol1").feature("t1").feature().remove("seDef");
    model.sol("sol1").attach("std1");
    model.sol("sol1").runAll();

    model.save('Resultfiles/Comsolmodel')
    return
