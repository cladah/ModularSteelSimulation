import logging

import numpy as np
from tc_python import *
from Framework.HelpFile import read_geninput
from Framework.Datastream_file import getaxisvalues


def rounding5(x, rnr=1):
    # Rounding to nearest 0.5
    base = .5/(10**rnr)
    return round(base * round(x/base), rnr+1)

def getTTTcompositions():
    """
    Returning a grid of compositions by varying carbon and nitrogen. Minimum step size of 0.5 w%

    :return: np.array of composition dicts
    """


    roundingTTT = 1
    ginput = read_geninput()
    TTTcompositions = list()
    fullcomposition = dict()
    for element in ginput['Material']['Composition'].keys():
        fullcomposition[element] = getaxisvalues("Composition/" + element)
    composition = ginput['Material']['Composition']

    mesh = list()
    # Creating grid of compositions
    for element in ginput['Material']['Composition'].keys():
        maxval = np.max(fullcomposition[element]) # changed from fullcomposition[element][-1]
        if composition[element] != maxval:
            if element == "C":
                tmplist = np.linspace(composition[element], maxval, 4)
                if rounding5(tmplist[0]) > 0:
                    tmplist = np.append([tmplist[0] - 0.05], tmplist)
                #if round(tmplist[0], 1) > 0:
                #    tmplist = np.concatenate([[tmplist[0] - 0.1], tmplist])
                tmplist = np.append(tmplist, [tmplist[-1] + 0.05])
                tmplist = [rounding5(elem) for elem in tmplist]
            elif element == "N":
                tmplist = np.linspace(composition[element], maxval, 2)
                if rounding5(tmplist[0]) > 0:
                    tmplist = np.append([tmplist[0] - 0.05], tmplist)
                tmplist = np.append(tmplist, [tmplist[-1] + 0.05])
                tmplist = [rounding5(elem) for elem in tmplist]
            else:
                tmplist = np.linspace(composition[element], maxval, 2)
                tmplist = [round(elem, 1) for elem in tmplist]

            tmplist = list(set(tmplist)) # getting unique values
            tmplist.sort()
        else:
            tmplist = [composition[element]]
        mesh.append(tmplist)

    g = np.meshgrid(*mesh)
    positions = np.vstack(list(map(np.ravel, g)))
    for compnr in range(len(positions[0, :])):  # The number 0 here is correlated to the coal as it varies the most
        tmpcomp = dict()
        i = 0
        for element in ginput['Material']['Composition'].keys():
            if element == "C" or element == "N":
                tmpcomp[element] = rounding5(positions[i, compnr])
            else:
                tmpcomp[element] = round(positions[i, compnr], 1)
            i = i + 1
        TTTcompositions.append(tmpcomp)
    return TTTcompositions

def TCequalibrium(ginput,minput,type):
    """

    :param type: Type of equalibrium

    :return: Activity for carbon and nitrogen, with graphite and gas as references.
    """
    print("Equalibrium calculation of type " + str(type))
    if type == "env":
        database = "SSUB6"
        material = setmaterial(ginput, minput, type)
        dependentmat = material[0]
        composition = material[1]
        phases = ["GAS", "C_S"]
        dormantphases = ["C_S"]
        referencestates = {"C": "C_S", "N": "GAS"} # N should have a different reference state? uses amonium?
    elif type == "mat":
        database = "TCFE12"
        dependentmat = ginput['Material']["Dependentmat"]
        composition = ginput['Material']["Composition"]
        phases = ["FCC_A1", "FCC_A1#2", "GAS", "GRAPHITE_A9"]
        dormantphases = ["GAS", "GRAPHITE_A9"]
        referencestates = {"C":"Graphite_A9", "N":"GAS"}
    else:
        raise KeyError

    with TCPython() as start:
        logging.getLogger("tc_python").setLevel(logging.ERROR)
        # create and configure a single equilibrium calculation
        calculation = (
            start
            .select_database_and_elements(database, [dependentmat] + list(composition))
            .get_system()
            .with_single_equilibrium_calculation()
            .set_condition(ThermodynamicQuantity.temperature(), minput["CNtemp"])
            .set_condition(ThermodynamicQuantity.pressure(), minput["CNPress"])
            .set_phase_to_suspended('*')
            #.disable_global_minimization()
        )

        for element in composition:
            calculation.set_condition(ThermodynamicQuantity.mass_fraction_of_a_component(element), composition[element]/100)
        for phase in phases:
            calculation.set_phase_to_entered(phase)
        for phase in dormantphases:
            calculation.set_phase_to_dormant(phase)
        for element in referencestates:
            calculation.with_reference_state(element, referencestates[element])
        logging.getLogger("tc_python").setLevel(logging.INFO)
        calc_result = (calculation.calculate())  # Aktiverar beräkningen
        activityC = calc_result.get_value_of(ThermodynamicQuantity.activity_of_component('C'))
        activityN = calc_result.get_value_of(ThermodynamicQuantity.activity_of_component('N'))
        print("Activity of carbon " + str(activityC))
        print("Activity of nitrogen " + str(activityN))
        if "C" in minput["DiffType"]:
            activityC = minput["Cpotential"]
        if "N" in minput["DiffType"]:
            activityN = activityN*500
            print("Changing nitrogen activity to " + str(activityN))
        return minput["Cpotential"], activityN

def TCDiffusionSolver(ginput, minput, Activity, compgrid):
    """
    :param activit: Activity of carbon and nitrogen [aC, aN]
    :param boosts: Nr of boost/rest cycles
    :param boost_t: Time with high activity at boundary
    :param diff_t: Time for diffusion between boosts
    :param compgrid: Composition grid for carbon and nitrogen
    :return: Coposition along x-axis.
    """
    print("Running Thermocalc diffusion solver")
    boosts = minput["BoostNr"]
    boost_t = minput["BoostTime"]
    diff_t = minput["DiffTime"]
    with TCPython() as session:
        logging.getLogger("tc_python").setLevel(logging.INFO)
        system = (session
                  .select_thermodynamic_and_kinetic_databases_with_elements("TCFE12", "MOBFE7",
                                                                            [ginput['Material']["Dependentmat"]] + list(ginput['Material']["Composition"]))
                  .without_default_phases().select_phase("FCC_A1").select_phase("GAS")
                  .select_phase("FCC_A1#2").select_phase("CEMENTITE_D011").select_phase("GRAPHITE_A9")
                  .get_system())
        """
        Adding information to the austenite profile.
        Material elements
        
        In order to read the composition for the last diffusion the composition must be read from the XDMF file.
        
        The last diffusion is for now overwritten.
        """
        profile0 = CompositionProfile()
        pos = getaxisvalues("nodes")
        for element in ginput['Material']["Composition"]:
            if element in ["C", "N"]:
                continue
            profile0.add(element, ElementProfile.constant(ginput['Material']["Composition"][element]))

        pointgrid = PointByPointGrid()
        for i in range(len(compgrid["C"])):
            tmppoint = GridPoint(pos[i][0])
            for element in ginput["Material"]["Composition"].keys():
                tmppoint.add_composition(element, compgrid[element][i])
            pointgrid.add_point(tmppoint)



        austenite = (Region("Austenite")
                     .with_point_by_point_grid_containing_compositions(pointgrid)
                     .add_phase("FCC_A1"))

        BC_Diffusion = BoundaryCondition.mixed_zero_flux_and_activity().set_zero_flux_for_element("C")
        #BC_Boost = BoundaryCondition.activity_flux_function().set_flux_function(element_name="C", f="-5E-8", n=1.0, g=str(Activity[1]))
        BC_Boost = BoundaryCondition.mixed_zero_flux_and_activity()
        if "N" in minput["DiffType"]:
            BC_Boost.set_activity_for_element("N", str(Activity[1]))
        if "C" in minput["DiffType"]:
            print("BOOST FOR C set to " + str(Activity[0]))
            BC_Boost = BoundaryCondition.activity_flux_function().set_flux_function(element_name="C", f="-5E-8", n=1.0,
                                                                                    g=str(Activity[0]))
            #BC_Boost.set_activity_for_element("C", str(Activity[0]))
        else:
            raise KeyError('Diffusion input wrong. Adjust The input to include Carbon (C) or Nitrogen (CN)')


        current_time = 0
        if len(minput["Temp"]) == 1:
            calc = system.with_isothermal_diffusion_calculation()
            calc.set_temperature(minput["Temp"][0])
        else:
            TempProf = TemperatureProfile()
            TempProf.add_time_temperature(current_time, minput["Temp"][0])
            for i in range(len(boost_t)):
                current_time = current_time + boost_t[i]
                TempProf.add_time_temperature(current_time, minput["Temp"][2*i])
                current_time = current_time + diff_t[i]
                TempProf.add_time_temperature(current_time, minput["Temp"][2*i+1])
            calc = system.with_non_isothermal_diffusion_calculation()
            calc.with_temperature_profile(TempProf)
        current_time = boost_t[0]
        #TemperatureProfile.add_time_temperature()
        #boost_calculation = (system.with_non_isothermal_diffusion_calculation().with_temperature_profile(TemperatureProfile.add_time_temperature()))
        boost_calculation = (calc
                             .with_reference_state("C", "GRAPHITE_A9")
                             .with_reference_state("N", "GAS")
                             .with_cylindrical_geometry()
                             .remove_all_regions()
                             .add_region(austenite)
                             .with_right_boundary_condition(BC_Boost)
                             .set_simulation_time(current_time))
        if current_time == 0:
            current_time = current_time + diff_t[0]
            diffusion_calculation = boost_calculation.with_right_boundary_condition(BC_Diffusion)
            diffusion_calculation.set_simulation_time(current_time)
        else:
            results_boost = boost_calculation.calculate()
            current_time = current_time + diff_t[0]
            # First diffusion
            diffusion_calculation = (results_boost.with_continued_calculation()
                                     .with_right_boundary_condition(BC_Diffusion)
                                     .set_simulation_time(current_time))
        results_diffusion = diffusion_calculation.calculate()
        for i in range(boosts):
            if i == 0:
                continue
            # Boost cycle
            current_time = current_time + boost_t[i]
            boost_calculation = (results_diffusion.with_continued_calculation()
                                 .with_right_boundary_condition(BC_Boost)
                                 .set_simulation_time(current_time))
            results_boost = boost_calculation.calculate()
            # Diffusion cycle
            current_time = current_time + diff_t[i]
            diffusion_calculation = (results_boost.with_continued_calculation()
                                     .with_right_boundary_condition(BC_Diffusion)
                                     .set_simulation_time(current_time))
            results_diffusion = diffusion_calculation.calculate()


        mass_frac = dict()
        composition = []
        composition.append(ginput['Material']["Dependentmat"])
        composition.append(ginput['Material']["Composition"])
        for element in ginput['Material']["Composition"]:
            distance, mass_frac_temp = results_diffusion.get_mass_fraction_of_component_at_time(element, SimulationTime.LAST)
            mass_frac[element] = mass_frac_temp
        return distance, mass_frac
def TCcarburizing_LPC(ginput, minput, activityair):
    print("Running Low pressure carburization model in ThermoCalc")
    """

    :param activityair: Activity of carbon and nitrogen [aC, aN]
    :param boosts: Nr of boost/rest cycles
    :param boost_t: Time with high activity at boundary
    :param rest_t: Time for diffusion between boosts
    :return: Coposition along x-axis.
    """
    data = read_input()
    boosts = minput["BoostNr"]
    boost_t = minput["BoostTime"]
    diff_t = minput["DiffTime"]
    with TCPython() as session:
        logging.getLogger("tc_python").setLevel(logging.INFO)
        system = (session
                  .select_thermodynamic_and_kinetic_databases_with_elements("TCFE12", "MOBFE7",
                                                                            [data['Material']["Dependentmat"]] + list(
                                                                                data['Material']["Composition"]))
                  .without_default_phases().select_phase("FCC_A1").select_phase("GAS").select_phase(
            "FCC_A1#2").select_phase("CEMENTITE_D011").select_phase("GRAPHITE_A9")
                  .get_system())
        """
        Adding information to the austenite profile.
        Material elements
        """

        austeniteprofile = CompositionProfile()
        for element in data['Material']["Composition"]:
            austeniteprofile.add(element, ElementProfile.constant(data['Material']["Composition"][element]))
            # Previous grid points PointByPointGrid.add_point()
            # Go through all the points for the previous simulation and add them to the grid.
            #
            # ElementProfile.funct(PointByPointGrid.add_point())
        austenite = (Region("Austenite").set_width(data['Geometry']["radius"])
                     .with_grid(CalculatedGrid.geometric()
                                .set_no_of_points(data['Geometry']["nodes"])
                                .set_geometrical_factor(data['Geometry']['meshscaling']))
                     .with_composition_profile(austeniteprofile)
                     .add_phase("FCC_A1"))

        BC_Diffusion = BoundaryCondition.mixed_zero_flux_and_activity().set_zero_flux_for_element("C")
        BC_Boost = BoundaryCondition.activity_flux_function().set_flux_function(element_name="C", f="-5E-8", n=1.0, g=str(1.0))

        current_time = boost_t[0]
        print(boost_t)
        boost_calculation = (system
                             .with_isothermal_diffusion_calculation()
                             .with_reference_state("C", "GRAPHITE_A9")
                             .set_temperature(data['Thermo']["CNtemp"])
                             .with_cylindrical_geometry()
                             .remove_all_regions()
                             .add_region(austenite)
                             .with_right_boundary_condition(BC_Boost)
                             .set_simulation_time(current_time))
        results_boost = boost_calculation.calculate()
        current_time = current_time + diff_t[0]
        # First diffusion
        diffusion_calculation = (results_boost.with_continued_calculation()
                                 .with_right_boundary_condition(BC_Diffusion)
                                 .set_simulation_time(current_time))
        results_diffusion = diffusion_calculation.calculate()
        for i in range(boosts):
            if i == 0:
                continue
            # Boost cycle
            current_time = current_time + boost_t[i]
            boost_calculation = (results_diffusion.with_continued_calculation()
                                 .with_right_boundary_condition(BC_Boost)
                                 .set_simulation_time(current_time))
            results_boost = boost_calculation.calculate()
            # Diffusion cycle
            current_time = current_time + diff_t[i]
            diffusion_calculation = (results_boost.with_continued_calculation()
                                     .with_right_boundary_condition(BC_Diffusion)
                                     .set_simulation_time(current_time))
            results_diffusion = diffusion_calculation.calculate()


        mass_frac = dict()
        composition = []
        composition.append(ginput["Dependentmat"])
        composition.append(ginput["Composition"])
        for element in ginput["Composition"]:
            distance, mass_frac_temp = results_diffusion.get_mass_fraction_of_component_at_time(element, SimulationTime.LAST)
            mass_frac[element] = mass_frac_temp
        return distance, mass_frac


def setmaterial(ginput,minput,type):
    if type == "mat":
        return ginput['Material']["Dependentmat"], ginput['Material']["Compositon"]
    elif type == "env":
        pass
    elif type == "Cenv":
        pass
    elif type == "Nenv":
        pass
    else:
        raise KeyError("TCcalculation error")

    if minput["CNenv"] == "Argon":
        print('Using argon as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 1, 'C': 1, 'O': 1}
    elif minput["CNenv"] == "Acetylene":
        """
        Low pressure atmosphere
        """
        print('Using acetylene as atmosphere')
        env_dep = 'N'
        #env_comp = {'H': 7.7, 'C': 92.2}
        env_comp = {'H': 3.7, 'C': 44.4}
    elif minput["CNenv"] == "Methane":
        print('Using methane as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 4.5788, 'C': 13.6344, 'O': 18.1614}
        activityair = [1.471, 0.639]
    elif minput["CNenv"] == 'Propane':
        print('Using propane as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 3.2762, 'C': 14.0238, 'O': 18.6801}
    else:
        raise KeyError("Wrong carbonitiding environment, use Argon, Methane, Acetylene, or Propane")

    return env_dep, env_comp

def calculatePearlite(temperatures, composition):
    print("Pearlite model")
    ginput = read_geninput()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = ginput['Material']["Dependentmat"]
    phases = ["FCC_A1"]

    with TCPython() as start:
        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kindatabase,
                                                                      [dependentmat] + list(composition))
            .deselect_phase("*")
            .select_phase("FCC_A1")
            .select_phase("BCC_A2")
            .select_phase("CEMENTITE")
            .get_system()
            .with_property_model_calculation("Pearlite").set_temperature(1000).set_composition_unit(
                CompositionUnit.MASS_PERCENT)
            # .set_argument()
        )
        for element in composition:
            calculation.set_composition(element, composition[element])

        # print("Available arguments: {}".format(calculation.get_arguments()))
        starttime, halftime, finishtime = list(), list(), list()
        calc_result = (calculation.set_temperature(temperatures[0])
                       .calculate()  # Aktiverar beräkningen
                       )
        for x in temperatures:
            calc_result = (calculation.set_temperature(x)
                           .calculate()  # Aktiverar beräkningen
                           )
            starttime.append(calc_result.get_value_of("Start"))
            halftime.append(calc_result.get_value_of("Half"))
            finishtime.append(calc_result.get_value_of("Finish"))
        print("Available result quantities: {}".format(calc_result.get_result_quantities()))
    return starttime, halftime, finishtime
        #'GrainSize', 'Criterion', 'PearliteMode', 'Austenite composition from', 'Austenitizing temperature'
        #for phase in phases:
        #    calculation.set_phase_to_entered(phase)
        #for phase in dormantphases:
        #    calculation.set_phase_to_dormant(phase)
        #for element in referencestates:
        #    calculation.with_reference_state(element,referencestates[element])

def calculateBainite(temperatures, composition):
    print("Bainite model")
    ginput = read_geninput()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = ginput['Material']["Dependentmat"]
    #composition = data['Material']["Composition"]
    phases = ["FCC_A1"]

    with TCPython() as start:
        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kindatabase,
                                                                      [dependentmat] + list(composition))
            .deselect_phase("*")
            .select_phase("FCC_A1")
            .select_phase("BCC_A2")
            .select_phase("CEMENTITE")
            .get_system()
            .with_property_model_calculation("Bainite").set_temperature(1000).set_composition_unit(CompositionUnit.MASS_PERCENT)
            # .set_argument()
        )
        for element in composition:
            calculation.set_composition(element, composition[element])

        #print("Available arguments: {}".format(calculation.get_arguments()))
        starttime, halftime, finishtime = list(), list(), list()
        for x in temperatures:
            calc_result = (calculation.set_temperature(x)
                           .calculate()  # Aktiverar beräkningen
                           )
            starttime.append(calc_result.get_value_of("Start time (2% bainite)"))
            halftime.append(calc_result.get_value_of("Half time (50% bainite)"))
            finishtime.append(calc_result.get_value_of("Finish time (98% bainite)"))
        #print("Available result quantities: {}".format(calc_result.get_result_quantities()))
    return starttime, halftime, finishtime
def calculateMartensite(composition):
    print("Martensite model")
    ginput = read_geninput()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = ginput['Material']["Dependentmat"]
    # composition = data['Material']["Composition"]
    phases = ["FCC_A1"]

    with TCPython() as start:
        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kindatabase,
                                                                      [dependentmat] + list(composition))
            .deselect_phase("*")
            .select_phase("FCC_A1")
            .select_phase("BCC_A2")
            .select_phase("CEMENTITE")
            .get_system()
            .with_property_model_calculation("Martensite Temperatures").set_temperature(1000).set_composition_unit(
                CompositionUnit.MASS_PERCENT)
            # .set_argument()
        )
        for element in composition:
            calculation.set_composition(element, composition[element])

        #print("Available arguments: {}".format(calculation.get_arguments()))
        calc_result = (calculation.calculate()  # Aktiverar beräkningen
                       )
        #print("Available result quantities: {}".format(calc_result.get_result_quantities()))
        start = calc_result.get_value_of("Ms")
        half = calc_result.get_value_of("M50")
        finish = calc_result.get_value_of("M99")
        #print([start, half, finish])
    return start, half, finish

def calculateFerrite(temperatures, composition):
    print("Ferrite model")
    ginput = read_geninput()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = ginput['Material']["Dependentmat"]
    #composition = data['Material']["Composition"]
    phases = ["FCC_A1"]
    tmpcomp = composition.copy()
    if tmpcomp["N"] == 0.0:
        print("Taking away nitrogen from calculation")
        tmpcomp.pop('N', None)
    with TCPython() as start:
        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kindatabase,
                                                                      [dependentmat] + list(tmpcomp))
            .deselect_phase("*")
            .select_phase("FCC_A1")
            .select_phase("BCC_A2")
            #.select_phase("CEMENTITE")
            .get_system()
            .with_property_model_calculation("Ferrite")
            .set_composition_unit(CompositionUnit.MASS_PERCENT)
            .set_argument("Growth mode", "Orthoequilibrium (OE)")
            .set_argument("GrainSize", "100")
            # .set_argument()
        )
        for element in tmpcomp:
            calculation.set_composition(element, tmpcomp[element])

        # print(calculation.get_model_description())
        # print("Available arguments: {}".format(calculation.get_arguments()))
        starttime, halftime, finishtime = list(), list(), list()
        for x in temperatures:
            try:
                calc_result = (calculation.set_temperature(x)
                               .calculate()  # Aktiverar beräkningen
                               )
                starttime.append(calc_result.get_value_of("Ferrite start"))
                halftime.append(calc_result.get_value_of("Ferrite half"))
                finishtime.append(calc_result.get_value_of("Ferrite finish"))
            except tc_python.exceptions.CalculationException:
                print("No result")
                starttime.append(np.nan)
                halftime.append(np.nan)
                finishtime.append(np.nan)
        # print("Available result quantities: {}".format(calc_result.get_result_quantities()))
    return starttime, halftime, finishtime
def calculateTTT(temperatures, composition):
    print("TTT model")
    ginput = read_geninput()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = ginput['Material']["Dependentmat"]
    # composition = data['Material']["Composition"]
    phases = ["FCC_A1"]

    with TCPython() as start:
        logging.getLogger("tc_python").setLevel(logging.ERROR)

        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kindatabase,
                                                                      [dependentmat] + list(composition))
            .deselect_phase("*")
            .select_phase("FCC_A1")
            .select_phase("BCC_A2")
            .select_phase("CEMENTITE")
            .get_system()
            .with_property_model_calculation("TTT Diagram").set_temperature(1000).set_composition_unit(
                CompositionUnit.MASS_PERCENT)
            # .set_argument()
        )
        # TTT Diagram - arguments
        # [2% austenite transformed, 50% austenite transformed, 98% austenite transformed, Bainite start (2%),
        # Ferrite start (2%), Martensite 50%, Martensite 98%, Martensite start,
        # Pearlite start (2%), Total ferrite start (2%), Total ferrite+cementite start (2%)]
        for element in composition:
            calculation.set_composition(element, composition[element])
        logging.getLogger("tc_python").setLevel(logging.INFO)
        calculation.set_argument("GrainSize", data["Thermo"]["GrainSize"])
        calculation.set_argument("Austenitizing temperature", data["Thermo"]["Austenitizingtemp"])
        print("Available arguments: {}".format(calculation.get_arguments()))
        aust2, aust50, aust98, mart_start, mart_half, \
            mart_end, pearlite, bainite, ferrite, ferrite_cem = list(), list(), list(), \
            list(), list(), list(), list(), list(), list(), list()
        legend = ["2% austenite transformed","50% austenite transformed","98% austenite transformed",
                  "Martensite start", "Martensite 50%", "Martensite 98%",
                  "Pearlite start (2%)", "Bainite start (2%)", "Ferrite start (2%)",
                  "Total ferrite+cementite start (2%)"]
        for x in temperatures:
            calc_result = (calculation.set_temperature(x)
                           .calculate()  # Aktiverar beräkningen
                           )
            aust2.append(calc_result.get_value_of("2% austenite transformed"))
            aust50.append(calc_result.get_value_of("50% austenite transformed"))
            aust98.append(calc_result.get_value_of("98% austenite transformed"))
            mart_start.append(calc_result.get_value_of("Martensite start"))
            mart_half.append(calc_result.get_value_of("Martensite 50%"))
            mart_end.append(calc_result.get_value_of("Martensite 98%"))
            pearlite.append(calc_result.get_value_of("Pearlite start (2%)"))
            bainite.append(calc_result.get_value_of("Bainite start (2%)"))
            ferrite.append(calc_result.get_value_of("Ferrite start (2%)"))
            ferrite_cem.append(calc_result.get_value_of("Total ferrite+cementite start (2%)"))
        TTT = [aust2, aust50, aust98, mart_start, mart_half, mart_end, pearlite, bainite, ferrite, ferrite_cem]
        print("Available result quantities: {}".format(calc_result.get_result_quantities()))
    return TTT, legend

def testTC():
    temperatures = [900, 800, 700, 600]
    print("Ferrite model")
    data = read_input()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    composition = data["Material"]["Composition"]
    dependentmat = data['Material']["Dependentmat"]
    # composition = data['Material']["Composition"]
    phases = ["FCC_A1"]

    with TCPython() as start:
        logging.getLogger("tc_python").setLevel(logging.ERROR)
        # [Martensite Temperatures, Interfacial energy, T-Zero Temperature, Crack Susceptibility Coefficient,
        # Pearlite, Spinodal, Critical Transformation Temperatures,
        # Bainite, Phase Transition, Equilibrium, CCT Diagram, Coarsening - Ni,
        # Ferrite, Equilibrium with Freeze-in Temperature,
        # Martensite Fractions, Liquidus and Solidus Temperature, Equilibrium with Freeze-in Temperature - Ni,
        # Antiphase Boundary Energy - Ni, Yield Strength, Driving Force, Coarsening, TTT Diagram]

        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kindatabase,
                                                                      [dependentmat] + list(composition))
            .deselect_phase("*")
            .select_phase("FCC_A1")
            .select_phase("BCC_A2")
            .select_phase("CEMENTITE")
            .get_system()
            .with_property_model_calculation("TTT Diagram").set_temperature(1000).set_composition_unit(
                CompositionUnit.MASS_PERCENT)
            # .set_argument()
        )
        # TTT Diagram - arguments
        # [2% austenite transformed, 50% austenite transformed, 98% austenite transformed, Bainite start (2%),
        # Ferrite start (2%), Martensite 50%, Martensite 98%, Martensite start,
        # Pearlite start (2%), Total ferrite start (2%), Total ferrite+cementite start (2%)]
        for element in composition:
            calculation.set_composition(element, composition[element])
        logging.getLogger("tc_python").setLevel(logging.INFO)
        calculation.set_argument("GrainSize", data["Thermo"]["GrainSize"])
        calculation.set_argument("Austenitizing temperature", data["Thermo"]["Austenitizingtemp"])
        print("Available arguments: {}".format(calculation.get_arguments()))
        starttime, halftime, finishtime = list(), list(), list()
        for x in temperatures:
            calc_result = (calculation.set_temperature(x)
                           .calculate()  # Aktiverar beräkningen
                           )
            starttime.append(calc_result.get_value_of("Ferrite start"))
            halftime.append(calc_result.get_value_of("Ferrite half"))
            finishtime.append(calc_result.get_value_of("Ferrite finish"))
        print("Available result quantities: {}".format(calc_result.get_result_quantities()))
    return starttime, halftime, finishtime
