import pathlib

import mph
import os
import csv
import numpy as np
from Sphere.HelpFile import readresultfile, read_input
from Sphere.Datastream_file import readdatastream, adjustdatastream, getaxisvalues
def modeldatatoComsolfiles():
    print("Adjusting phase transformation data to Comsol specifics")
    xyz = readdatastream("nodes")
    phases = ["Ferrite", "Perlite", "Bainite", "Martensite"]
    for phase in phases:
        if phase != "Martensite":
            tmpT = readresultfile("Modeldata", phase + "/JMAK/T")
            tmp1 = readresultfile("Modeldata", phase + "/JMAK/tau")
            tmp2 = readresultfile("Modeldata", phase + "/JMAK/n")
            tmpnames = ["tau", "n"]
        else:
            tmp1 = readresultfile("Modeldata", phase + "/KM/Ms")
            tmp2 = readresultfile("Modeldata", phase + "/KM/beta")
            tmpnames = ["Ms", "beta"]

        for i in [0, 1]:
            with open("Resultfiles/" + phase + "_" + tmpnames[i] + ".csv", 'w') as file:
                writer = csv.writer(file)
                r = [np.sqrt(xyz[j][0]**2 + xyz[j][1]**2) for j in range(len(xyz) - 1)]
                if phase in ["Ferrite", "Perlite", "Bainite"]:
                    if not tmp1.all():
                        break
                    if i == 0:
                        for j in range(len(tmp1) - 1):
                            for k in range(len(tmpT) - 1):
                                writer.writerow([r[j], tmpT[k], tmp1[j][k]])
                    else:
                        for j in range(len(tmp2) - 1):
                            for k in range(len(tmpT) - 1):
                                writer.writerow([r[j], tmpT[k], tmp2[j][k]])
                else:

                    if i == 0:
                        for j in range(len(tmp1) - 1):
                            writer.writerow([r[j], tmp1[j][0]])
                    else:
                        for j in range(len(tmp2) - 1):
                            writer.writerow([r[j], tmp2[j][0]])
        print(str(phase) + " data written to csv files")
    pass
def setupComsolSolver(model):
    print("Setting upp Comsol model")
    model.sol("sol1").study("std1")
    model.study("std1").feature("time").set("notlistsolnum", 1)
    model.study("std1").feature("time").set("notsolnum", "auto")
    model.study("std1").feature("time").set("listsolnum", 1)
    model.study("std1").feature("time").set("solnum", "auto")
    model.sol("sol1").feature().remove("t1")
    model.sol("sol1").feature().remove("v1")
    model.sol("sol1").feature().remove("st1")
    model.sol("sol1").create("st1", "StudyStep")
    model.sol("sol1").feature("st1").set("study", "std1")
    model.sol("sol1").feature("st1").set("studystep", "time")
    model.sol("sol1").create("v1", "Variables")
    model.sol("sol1").feature("v1").feature("comp1_audc_phase1_xi").set("scalemethod", "manual")
    model.sol("sol1").feature("v1").feature("comp1_audc_etripc13").set("scalemethod", "manual")
    model.sol("sol1").feature("v1").feature("comp1_audc_etripc11").set("scalemethod", "manual")
    model.sol("sol1").feature("v1").feature("comp1_audc_etripc33").set("scalemethod", "manual")
    model.sol("sol1").feature("v1").feature("comp1_audc_phase5_xi").set("scalemethod", "manual")
    model.sol("sol1").feature("v1").feature("comp1_u").set("scalemethod", "manual")
    model.sol("sol1").feature("v1").feature("comp1_audc_phase1_xi").set("scaleval", "1")
    model.sol("sol1").feature("v1").feature("comp1_audc_etripc13").set("scaleval", "1e-3")
    model.sol("sol1").feature("v1").feature("comp1_audc_etripc11").set("scaleval", "1e-3")
    model.sol("sol1").feature("v1").feature("comp1_audc_etripc33").set("scaleval", "1e-3")
    model.sol("sol1").feature("v1").feature("comp1_audc_phase5_xi").set("scaleval", "1")
    model.sol("sol1").feature("v1").feature("comp1_u").set("scaleval", "1e-2*0.01414213562373095")
    model.sol("sol1").feature("v1").set("control", "time")
    model.sol("sol1").create("t1", "Time")
    model.sol("sol1").feature("t1").set("tlist", "range(0,0.1,10)")
    model.sol("sol1").feature("t1").set("plot", "off");
    model.sol("sol1").feature("t1").set("plotgroup", "pg1");
    model.sol("sol1").feature("t1").set("plotfreq", "tout");
    model.sol("sol1").feature("t1").set("probesel", "all");
    model.sol("sol1").feature("t1").set("probes", ""); # -----------------------------------
    model.sol("sol1").feature("t1").set("probefreq", "tsteps");
    model.sol("sol1").feature("t1").set("rtol", 0.001);
    model.sol("sol1").feature("t1").set("atolglobalvaluemethod", "factor");
    model.sol("sol1").feature("t1").set("atolmethod", "comp1_audc_etripc11", "global", "comp1_audc_etripc13", "global", "comp1_audc_etripc33", "global", "comp1_audc_phase1_xi", "global", "comp1_audc_phase5_xi", "global", "comp1_T", "global", "comp1_u", "global");
    model.sol("sol1").feature("t1").set("atol", "comp1_audc_etripc11", "1e-3", "comp1_audc_etripc13", "1e-3", "comp1_audc_etripc33", "1e-3", "comp1_audc_phase1_xi", "1e-3", "comp1_audc_phase5_xi", "1e-3", "comp1_T", "1e-3", "comp1_u", "1e-3");
    model.sol("sol1").feature("t1").set("atolvaluemethod", "comp1_audc_etripc11", "factor", "comp1_audc_etripc13", "factor", "comp1_audc_etripc33", "factor", "comp1_audc_phase1_xi", "factor", "comp1_audc_phase5_xi", "factor", "comp1_T", "factor", "comp1_u", "factor");
    model.sol("sol1").feature("t1").set("atolvaluemethod", "comp1_audc_etripc11", "factor", "comp1_audc_etripc13", "factor", "comp1_audc_etripc33", "factor","comp1_audc_phase1_xi", "factor", "comp1_audc_phase5_xi", "factor", "comp1_T", "factor", "comp1_u", "factor");
    model.sol("sol1").feature("t1").set("reacf", True);
    model.sol("sol1").feature("t1").set("storeudot", True);
    model.sol("sol1").feature("t1").set("endtimeinterpolation", True);
    model.sol("sol1").feature("t1").set("estrat", "exclude");
    model.sol("sol1").feature("t1").set("maxorder", 2);
    model.sol("sol1").feature("t1").set("control", "time");
    model.sol("sol1").feature("t1").feature("aDef").set("cachepattern", True);
    model.sol("sol1").feature("t1").create("seDef", "Segregated");
    model.sol("sol1").feature("t1").create("se1", "Segregated");
    model.sol("sol1").feature("t1").feature("se1").feature().remove("ssDef");
    model.sol("sol1").feature("t1").feature("se1").create("ss1", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("segvar", "comp1_T", "comp1_u");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("subdamp", 1);
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("subjtech", "once");
    model.sol("sol1").feature("t1").create("d1", "Direct");
    model.sol("sol1").feature("t1").feature("d1").label("Direct (merged)");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").set("linsolver", "d1");
    model.sol("sol1").feature("t1").feature("se1").feature("ss1").label("Merged variables");
    model.sol("sol1").feature("t1").feature("se1").create("ss2", "SegregatedStep");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("segvar", "comp1_audc_etripc11", "comp1_audc_etripc13", "comp1_audc_etripc33", "comp1_audc_phase1_xi","comp1_audc_phase5_xi");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").set("linsolver", "d1");
    model.sol("sol1").feature("t1").feature("se1").feature("ss2").label("Austenite Decomposition");
    model.sol("sol1").feature("t1").feature("se1").set("segstabacc", "segaacc");
    model.sol("sol1").feature("t1").feature("se1").set("segaaccdim", 5);
    model.sol("sol1").feature("t1").feature("se1").set("segaaccdelay", 0);
    model.sol("sol1").feature("t1").feature("se1").set("segaaccmix", 0.9);
    model.sol("sol1").feature("t1").feature("se1").create("ll1", "LowerLimit");
    model.sol("sol1").feature("t1").feature("se1").feature("ll1").set("lowerlimit", "root.comp1.audc.phase5.xi 0 root.comp1.audc.phase1.xi 0 comp1.T 0 ");


    model.sol("sol1").feature("t1").feature("se1").create("ul1", "UpperLimit");
    model.sol("sol1").feature("t1").feature("se1").feature("ul1").set("upperlimit", " root.comp1.audc.phase1.xi 1 root.comp1.audc.phase5.xi 1");
    model.sol("sol1").feature("t1").create("i1", "Iterative");
    model.sol("sol1").feature("t1").feature("i1").set("linsolver", "gmres");
    model.sol("sol1").feature("t1").feature("i1").set("prefuntype", "left");
    model.sol("sol1").feature("t1").feature("i1").set("itrestart", 50);
    model.sol("sol1").feature("t1").feature("i1").set("rhob", 20);
    model.sol("sol1").feature("t1").feature("i1").set("maxlinit", 10000);
    model.sol("sol1").feature("t1").feature("i1").set("nlinnormuse", "on");
    model.sol("sol1").feature("t1").feature("i1").label("AMG, heat transfer variables ht (te1)");
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
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("i1").feature("mg1").feature("cs").feature("d1").set("pivotperturb", 1.0E-13);
    model.sol("sol1").feature("t1").feature().remove("fcDef");
    model.sol("sol1").feature("t1").feature().remove("seDef");
    return model
def setupComsol(model):
    directory = os.getcwd()
    if os.name == "posix":
        meshdirec = directory + '/Resultfiles/Mesh.nas'
        savedirec = directory + '/Resultfiles'
    else:
        meshdirec = directory + '/Resultfiles/Mesh.nas'
        savedirec = directory + '/Resultfiles'

    model.modelNode().create("comp1", True)

    # --------------- Setting up geometry and mesh ------------------#
    model.component("comp1").geom().create("geom1", 2)
    #model.component("comp1").geom("geom1").axisymmetric(True)
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
        model.func().create("h_" + ph, "Interpolation")
        model.nodeGroup(phase).add("func", "h_" + ph)
        model.func().create("alpha_k_"+ph, "Interpolation")
        model.nodeGroup(phase).add("func", "alpha_k_"+ph)
        if phase == "Austenite":
            pass
        elif phase == "Martensite":
            model.func().create("Ms_"+ph, "Interpolation")
            model.nodeGroup(phase).add("func", "Ms_"+ph)
            model.func().create("beta_"+ph, "Interpolation")
            model.nodeGroup(phase).add("func", "beta_"+ph)
        else:
            model.func().create("tau_"+ph, "Interpolation")
            model.nodeGroup(phase).add("func", "tau_"+ph)
            model.func().create("n_"+ph, "Interpolation")
            model.nodeGroup(phase).add("func", "n_"+ph)

    model.func().create("htc", "Interpolation")
    # --------------- Setting up material ------------------#
    for phase in phases:
        ph = phase[0]
        model.component("comp1").material().create(phase, "Common")
        model.component("comp1").material(phase).propertyGroup("def").set("heatcapacity", "Cp_" + ph + "(T)")
        model.component("comp1").material(phase).propertyGroup("def").set("thermalconductivity", "k_" + ph + "(T)")
        model.component("comp1").material(phase).propertyGroup().create("ThermalExpansion", "Thermal_expansion")
        model.component("comp1").material(phase).propertyGroup("ThermalExpansion").set("thermalexpansioncoefficient", "alpha_k_" + ph + "(T)")
        model.component("comp1").material(phase).propertyGroup().create("Enu", "Young's_modulus_and_Poisson's_ratio")
        model.component("comp1").material(phase).propertyGroup("Enu").set("E", "E_" + ph + "(T)")
        model.component("comp1").material(phase).propertyGroup("Enu").set("density", "7800[kg/m^3]")
        model.component("comp1").material(phase).propertyGroup("Enu").set("nu", "0.3")
        model.component("comp1").material(phase).propertyGroup().create("ElastoplasticModel", "Elastoplastic material model")
        model.component("comp1").material(phase).propertyGroup("ElastoplasticModel").set("sigmags", "Sy_" + ph + "(T)")

    # --------------- Adding physics ------------------#
    model.component("comp1").physics().create("solid", "SolidMechanics", "geom1")
    model.component("comp1").physics().create("ht", "HeatTransfer", "geom1")
    model.component("comp1").physics().create("audc", "AusteniteDecomposition", "geom1")
    model.component("comp1").physics("solid").prop("StructuralTransientBehavior").set("StructuralTransientBehavior",
                                                                                      "Quasistatic")

    model.component("comp1").physics("audc").prop("MaterialProperties").runCommand("makecompoundmaterial")
    model.component("comp1").material("audcmat").selection().all()
    model.component("comp1").material("audcmat").propertyGroup("def").set("thermalexpansioncoefficient", 22.0E-6)

    # Austenite
    model.component("comp1").physics("audc").feature("phase1").set("phaseMaterial", "Austenite")

    # Ferrite
    model.component("comp1").physics("audc").feature("phase2").set("phaseMaterial", "Ferrite")
    model.component("comp1").physics("audc").feature("phase2").selection().all()
    model.component("comp1").physics("audc").feature("ptran1").set("ptModel", "JMAK")
    model.component("comp1").physics("audc").feature("ptran1").set("taujmak", "tau_F(sqrt(x^2+y^2),T)")
    model.component("comp1").physics("audc").feature("ptran1").set("njmak", "n_F(sqrt(x^2+y^2),T)")
    #model.component("comp1").physics("audc").feature("ptran1").set("trip", True)

    # Perlite
    model.component("comp1").physics("audc").feature("phase3").set("phaseMaterial", "Perlite")
    model.component("comp1").physics("audc").feature("phase3").selection().all()
    model.component("comp1").physics("audc").feature("ptran2").set("ptModel", "JMAK")
    model.component("comp1").physics("audc").feature("ptran2").set("taujmak", "tau_P(sqrt(x^2+y^2),T)")
    model.component("comp1").physics("audc").feature("ptran2").set("njmak", "n_P(sqrt(x^2+y^2),T)")
    #model.component("comp1").physics("audc").feature("ptran2").set("trip", True)

    # Bainite
    model.component("comp1").physics("audc").feature("phase4").set("phaseMaterial", "Bainite")
    model.component("comp1").physics("audc").feature("phase4").selection().all()
    model.component("comp1").physics("audc").feature("ptran3").set("ptModel", "JMAK")
    model.component("comp1").physics("audc").feature("ptran3").set("taujmak", "tau_B(sqrt(x^2+y^2),T)")
    model.component("comp1").physics("audc").feature("ptran3").set("njmak", "n_B(sqrt(x^2+y^2),T)")
    #model.component("comp1").physics("audc").feature("ptran3").set("trip", True)

    # Martensite
    model.component("comp1").physics("audc").feature("phase5").set("phaseMaterial", "Martensite")
    model.component("comp1").physics("audc").feature("phase5").selection().all()
    model.component("comp1").physics("audc").feature("ptran4").set("Ms", "Ms_M(sqrt(x^2+y^2))")
    model.component("comp1").physics("audc").feature("ptran4").set("beta", "beta_M(sqrt(x^2+y^2))")

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
    model.component("comp1").physics("solid").create("symp1", "SymmetrySolid", 1)
    model.component("comp1").physics("solid").feature("symp1").selection().set(1, 2)
    model.component("comp1").physics("solid").feature("symp1").set("constraintMethod", "nodal")
    model.component("comp1").physics("ht").create("symp1", "Symmetry", 1)
    model.component("comp1").physics("ht").feature("symp1").selection().set(1, 2)
    #model.component("comp1").physics("ht").create("temp1", "TemperatureBoundary", 1)
    #model.component("comp1").physics("ht").feature("temp1").selection().set(3)
    #model.component("comp1").physics("ht").feature("temp1").set("T0", "100[degC]")
    # Boundary heat flux
    model.component("comp1").physics("ht").create("hf1", "HeatFluxBoundary", 1)
    model.component("comp1").physics("ht").feature("hf1").set("HeatFluxType", "ConvectiveHeatFlux")
    model.component("comp1").physics("ht").feature("hf1").selection().set(3)
    model.component("comp1").physics("ht").feature("hf1").set("h", "htc(T)")
    model.component("comp1").physics("ht").feature("hf1").set("Text", "30 [degC]")
    #model.component("comp1").physics("ht").create("hf1", "HeatFluxBoundary", 1)
    #model.component("comp1").physics("ht").feature("hf1").selection().set(3)
    #model.component("comp1").physics("ht").feature("hf1").set("q0_input", -1234.)

    #model.component("comp1").physics("ht").feature("hf1").set("q0_input", 12345)  # Add the correct temperature inflow

    #model.component("comp1").physics("ht").create("sym1", "Symmetry", 1)
    #model.component("comp1").physics("ht").feature("sym1").selection().set(2)




    # --------------- Setting up phase transformation ------------------#
    model.study().create("std1");
    model.study("std1").create("time", "Transient")
    model.study("std1").feature("time").setSolveFor("/physics/solid", True)
    model.study("std1").feature("time").setSolveFor("/physics/ht", True)
    model.study("std1").feature("time").setSolveFor("/physics/audc", True)
    model.study("std1").feature("time").setSolveFor("/multiphysics/lht1", True)
    model.study("std1").feature("time").setSolveFor("/multiphysics/ptstr1", True)
    model.study("std1").feature("time").setSolveFor("/multiphysics/te1", True)
    # model.func().create("an1", "Analytic")
    # model.func("an1").set("funcname", "Ms")
    # model.func("an1").set("expr", "0.011+0.01*carb[1/K]")
    # model.func("an1").set("args", "carb")
    # model.func("an1").set("expr", "200+100*carb[degC]")

    # --------------- Creating study ------------------#
    #model.component("comp1").physics("audc").active(False)
    #model.component("comp1").multiphysics("lht1").active(False)
    #model.component("comp1").multiphysics("ptstr1").active(False)



    #model = setupComsolSolver(model)
    model.save('Resultfiles/Comsolmodel')
    return model

def adjustComsol(model):
    print("Importing mesh to Comsol model")
    model.component("comp1").mesh("mesh1").feature("imp1").importData()
    model.component("comp1").mesh("mesh1").run()

    print("Adding phase transformation models to Comsol model")
    #JMAK_B = readresultfile("TTT_surface.hdf5","Bainite/JMAK")
    #JMAK_P = readresultfile("TTT_surface.hdf5","Perlite/JMAK")
    data = read_input()

    # Heat flux
    x = data["FEM"]["heatflux"]["T"]
    y = data["FEM"]["heatflux"]["htc"]
    for i in range(len(x)):
        model.func("htc").setIndex("table", x[i], i, 0)
        model.func("htc").setIndex("table", y[i], i, 1)
    model.func("htc").setIndex("argunit", "degC", 0)
    model.func("htc").setIndex("fununit", "W/(m^2*K)", 0)

    # Carbonnitriding temperature
    model.component("comp1").physics("ht").feature("init1").set("Tinit", "900[degC]")
    materialprop = ["E","Cp","k","Sy", "alpha_k", "h"]
    materials = ["Austenite","Ferrite","Perlite","Bainite","Martensite"]
    propunit = ["GPa", "J/(kg*K)", "W/(m*K)", "MPa", "1/K", "GPa"] # h = GPa
    for mat in materials:
        if mat in ["Ferrite","Perlite","Bainite"]:
            materialprop = ["E", "Cp", "k", "Sy", "alpha_k", "h", "tau", "n"]
        elif mat == "Martensite":
            materialprop = ["E", "Cp", "k", "Sy", "alpha_k", "h", "Ms","beta"]
        for prop in materialprop:
            # model.func("E_A").set("table", "[]")
            if prop in ["tau","n"]:
                model.func(prop + "_" + mat[0]).set("source", "file")
                model.func(prop + "_" + mat[0]).set("filename", "Resultfiles/" + mat + "_" + prop + ".csv")
                model.func(prop + "_" + mat[0]).set("nargs", "2")
                model.func(prop + "_" + mat[0]).setIndex("argunit", "m", 0)
                model.func(prop + "_" + mat[0]).setIndex("argunit", "K", 1)
            elif prop == "Ms":
                model.func(prop + "_" + mat[0]).set("source", "file")
                model.func(prop + "_" + mat[0]).set("filename", "Resultfiles/" + mat + "_" + prop + ".csv")
                model.func(prop + "_" + mat[0]).set("nargs", "1")
                model.func(prop + "_" + mat[0]).setIndex("argunit", "m", 0)
                model.func(prop + "_" + mat[0]).setIndex("fununit", "K", 1)
            elif prop == "beta":
                model.func(prop + "_" + mat[0]).set("source", "file")
                model.func(prop + "_" + mat[0]).set("filename", "Resultfiles/" + mat + "_" + prop + ".csv")
                model.func(prop + "_" + mat[0]).set("nargs", "1")
                model.func(prop + "_" + mat[0]).setIndex("argunit", "m", 0)
                model.func(prop + "_" + mat[0]).setIndex("fununit", "1/K", 0)

            else:
                x = data["Material"][mat][prop]["T"]
                y = data["Material"][mat][prop][prop]
                for i in range(len(x)):
                    model.func(prop + "_" + mat[0]).setIndex("table", x[i], i, 0)
                    model.func(prop + "_" + mat[0]).setIndex("table", y[i], i, 1)
                model.func(prop + "_" + mat[0]).setIndex("argunit", "degC", 0)
                model.func(prop + "_" + mat[0]).setIndex("fununit", propunit[materialprop.index(prop)], 0)


    model.save('Resultfiles/Comsolmodel')
    return model

def Comsolexport(model):

    resultdata = ["solid.eel11", "solid.eel12", "solid.eel22", "solid.eel23", "solid.eel13", "solid.eel33",
                  "solid.sl11", "solid.sl12", "solid.sl22", "solid.sl23", "solid.sl13", "solid.sl33",
                  "T", "audc.phase1.xi", "audc.phase2.xi", "audc.phase3.xi", "audc.phase4.xi", "audc.phase5.xi"]
    resultnames = ["eel11", "eel12", "eel22", "eel23", "eel13", "eel33", "sl11", "sl12", "sl22", "sl23", "sl13",
                   "sl33", "T", "Austenite", "Ferrite", "Perlite", "Bainite", "Martensite"]
    model.result().export().create("data1", "Data")
    model.result().export("data1").set("filename", "tmpComsol.csv")
    model.result().export("data1").setIndex("looplevelinput", "all", 0)
    # model.result().export("data1").set("header", False)
    data_dict = {}

    for i in range(len(resultdata)):
        model.result().export("data1").setIndex("expr", resultdata[i], 0)
        model.result().export("data1").run()


        with open(pathlib.Path(".").absolute() / "Resultfiles" / "tmpComsol.csv") as file:
            csvreader = csv.reader(file)
            for j in range(8):
                next(csvreader)

            time = next(csvreader)[2:]
            time = [j.split("=")[1] for j in time]
            j = 0
            data = []
            x = []
            y = []

            # Get index
            for line in csvreader:
                x.append(line[0])
                y.append(line[1])
                data.append(line[2:])
                j = j + 1
            data = np.array(data).astype(float)
            indx = getComsolindx(x, y)


        for j in range(len(time)):
            if str(time[j]) not in data_dict.keys():
                data_dict[str(time[j])] = dict()
            data_dict[str(time[j])][resultnames[i]] = data[indx, j]

    for i in range(len(time)):
        adjustdatastream(data_dict[str(time[i])], t_data=time[i])
        print("Exported timestep " + time[i])
def getComsolindx(xdata, ydata):
    xdata = np.array(xdata).astype(float)
    ydata = np.array(ydata).astype(float)
    nodes = readdatastream("nodes")
    xdata = np.around(xdata, 8)
    ydata = np.around(ydata, 8)
    x = np.around(np.array(nodes[:, 0]), 8)
    y = np.around(np.array(nodes[:, 1]), 8)

    # Length correct?

    # Getting position of comsol nodes
    indxcomsol = np.ones(len(nodes))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if x[j] == xdata[i]:
                if y[j] == ydata[i]:
                    indxcomsol[j] = i
                    continue
    indxcomsol = indxcomsol.astype(int)
    return indxcomsol

def addComsoldatadatastream():
    directory = os.getcwd()
    cachedirectory = directory + '/Cachefiles'
    savedirec = directory + '/Resultfiles'
    data = np.loadtxt(savedirec + "/Phasecomp.txt")
    nodes = readdatastream("nodes")

    xdata = np.around(np.array(data)[:, 0], 8)
    ydata = np.around(np.array(data)[:, 1], 8)
    x = np.around(np.array(nodes[:, 0]), 8)
    y = np.around(np.array(nodes[:, 1]), 8)

    # Getting position of comsol nodes
    indxcomsol = np.ones(len(nodes))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if x[j] == xdata[i]:
                if y[j] == ydata[i]:
                    indxcomsol[j] = i
                    continue
    indxcomsol = indxcomsol.astype(int)
    t = [*range(0, 30, 1), *range(60, 600, 60)]
    adjustdatastream("Austenite", data[:, 2][indxcomsol], "nodes")


def resultconverter():
    import numpy as np
    directory = os.getcwd()
    cachedirectory = directory + '/Cachefiles'
    savedirec = directory + '/Resultfiles'
    nodes = readdatastream("nodes")

    data = np.loadtxt(savedirec + "/Phasecomp.txt")
    xdata = np.around(np.array(data)[:, 0], 8)
    ydata = np.around(np.array(data)[:, 1], 8)
    x = np.around(np.array(nodes[:, 0]), 8)
    y = np.around(np.array(nodes[:, 1]), 8)

    # Getting position of comsol nodes
    indxcomsol = np.ones(len(nodes))
    for i in range(len(nodes)):
        for j in range(len(nodes)):
            if x[j] == xdata[i]:
                if y[j] == ydata[i]:
                    indxcomsol[j] = i
                    continue
    indxcomsol = indxcomsol.astype(int)

    adjustdatastream("Austenite", data[:, 2][indxcomsol], "nodes")
    adjustdatastream("Ferrite", data[:, 3][indxcomsol], "nodes")
    adjustdatastream("Pearlite", data[:, 4][indxcomsol], "nodes")
    adjustdatastream("Bainite", data[:, 5][indxcomsol], "nodes")
    adjustdatastream("Martensite", data[:, 6][indxcomsol], "nodes")

    data = np.loadtxt(savedirec + "/Stress.txt")
    adjustdatastream("vonMises", data[:, 2][indxcomsol], "nodes")

    data = np.loadtxt(savedirec + "/Strain.txt")
    adjustdatastream("ep1", data[:, 2][indxcomsol], "nodes")
    adjustdatastream("ep2", data[:, 3][indxcomsol], "nodes")
    adjustdatastream("ep3", data[:, 4][indxcomsol], "nodes")
    adjustdatastream("evol", data[:, 5][indxcomsol], "nodes")
    adjustdatastream("edeve", data[:, 6][indxcomsol], "nodes")

def runComsol(parent):
    directory = os.getcwd()
    savedirec = directory + '/Resultfiles'
    client = mph.start()
    #pymodel = client.load("Resultfiles/Comsolmodel.mph")
    #model = pymodel.java
    modeldatatoComsolfiles()
    parent.updateprogress(0.1)
    print("Opening Comsol multiphysics")
    #client = mph.start()

    print("Setting up model")
    if os.path.isfile("Resultfiles/Comsolmodel.mph"):
        print("Comsolmodel.mph file exist")
        pymodel = client.load("Resultfiles/Comsolmodel.mph")
        model = pymodel.java
    else:
        print("Comsolmodel.mph file doesn't exist")
        print("Creating Comsol model")
        pymodel = client.create()
        model = pymodel.java
        model = setupComsol(model)
    #model.util.ModelUtil.showProgress(savedirec + "/Comsolprogress.txt")
    print("Adjusting model to input")
    model = adjustComsol(model)
    parent.updateprogress(0.3)
    print("Running model")
    model.study("std1").feature("time").set("tlist", "range(0,1,30),range(60,60,600)")
    model.study("std1").run()
    model.save('Resultfiles/Comsolmodel')
    parent.updateprogress(0.9)
    print("Comsol model successfully ran")
    print("Exporting results")
    Comsolexport(model)

    client.clear()
    #resultconverter()
    parent.updateprogress(1.0)
    print("Comsol multiphysics closed")


    #client = mph.start()
    #pymodel = client.load("Resultfiles/Comsolmodel.mph")
    #model = pymodel.java

    #model.save('Resultfiles/Comsolmodel')

    #setupComsolSolver()