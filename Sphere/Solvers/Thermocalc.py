import numpy as np
from tc_python import *

def TCequalibrium(type):
    from HelpFile import read_input
    data = read_input()
    if type == "env":
        database = "SSUB6"
        material = setmaterial(data, type)
        dependentmat = material[0]
        composition = material[1]
        phases = ["GAS", "C_S"]
        dormantphases = ["C_S"]
        referencestates = {"C": "C_S", "N": "GAS"}
    elif type == "mat":
        database = "TCFE12"
        dependentmat = data['Material']["Dependentmat"]
        composition = data['Material']["Composition"]
        phases = ["FCC_A1", "FCC_A1#2", "GAS", "GRAPHITE_A9"]
        dormantphases = ["GAS", "GRAPHITE_A9"]
        referencestates = {"C":"Graphite_A9", "N":"GAS"}
    else:
        raise KeyError

    with TCPython() as start:
        # create and configure a single equilibrium calculation
        calculation = (
            start
            .select_database_and_elements(database, [dependentmat] + list(composition))
            .get_system()
            .with_single_equilibrium_calculation()
            .set_condition(ThermodynamicQuantity.temperature(), data['Thermo']["CNtemp"])
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
            calculation.with_reference_state(element,referencestates[element])
        calc_result = (calculation
                       .calculate()  # Aktiverar beräkningen
                       )
        activityC = calc_result.get_value_of(ThermodynamicQuantity.activity_of_component('C'))
        activityN = calc_result.get_value_of(ThermodynamicQuantity.activity_of_component('N'))
        return activityC, activityN

def TCcarbonitriding(activityair):
    from HelpFile import read_input
    data = read_input()
    with TCPython() as session:
        system = (session
                  .select_thermodynamic_and_kinetic_databases_with_elements("TCFE12", "MOBFE7", [data['Material']["Dependentmat"]] + list(data['Material']["Composition"]))
                  .without_default_phases().select_phase("FCC_A1").select_phase("GAS").select_phase("FCC_A1#2").select_phase("CEMENTITE_D011").select_phase("GRAPHITE_A9")
                  .get_system())
        austenite = Region("Austenite")
        austenite.set_width(data['Geometry']["radius"])
        austenite.with_grid(CalculatedGrid.geometric()
                       .set_no_of_points(data['Geometry']["nodes"])
                       .set_geometrical_factor(data['Geometry']['meshscaling']))
        austenite.add_phase("FCC_A1")
        austenite.add_phase("FCC_A1#2")
        austeniteprofile = CompositionProfile()
        for element in data['Material']["Composition"]:
            austeniteprofile.add(element, ElementProfile.constant(data['Material']["Composition"][element]))
        austenite.with_composition_profile(austeniteprofile)
        calculation = (system
                      .with_isothermal_diffusion_calculation()
                      .with_reference_state("N", "GAS").with_reference_state("C", "GRAPHITE_A9")
                      .set_temperature(data['Thermo']["CNtemp"])
                      .set_simulation_time(data['Thermo']["CNtime"] * 3600)
                      .with_right_boundary_condition(BoundaryCondition.mixed_zero_flux_and_activity()
                                                     .set_activity_for_element('C', str(activityair[0]))
                                                     .set_activity_for_element('N', str(activityair[1])))
                      .with_spherical_geometry().remove_all_regions()
                      .add_region(austenite))


        result = calculation.calculate()
        mass_frac = dict()
        composition = []
        composition.append(data['Material']["Dependentmat"])
        composition.append(data['Material']["Composition"])
        for element in data['Material']["Composition"]:
            distance = []
            distance, mass_frac_temp = result.get_mass_fraction_of_component_at_time(element, SimulationTime.LAST)
            mass_frac[element] = mass_frac_temp
        return distance, mass_frac

def setmaterial(data,type):
    if type == "mat":
        return data['Material']["Dependentmat"],data['Material']["Compositon"]
    elif type == "env":
        pass
    else:
        raise KeyError("TCcalculation error")

    if data['Thermo']["CNenv"] =="Argon":
        print('Using argon as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 1, 'C': 1, 'O': 1}
    elif data['Thermo']["CNenv"] == "Methane":
        print('Using methane as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 4.5788, 'C': 13.6344, 'O': 18.1614}
        activityair = [1.471, 0.639]
    elif data['Thermo']["CNenv"] =='Propane':
        print('Using propane as atmosphere')
        env_dep = 'N'
        env_comp = {'H': 3.2762, 'C': 14.0238, 'O': 18.6801}
    else:
        raise KeyError("Wrong carbonitiding environment, use Argon, Methane, or Propane")

    return env_dep, env_comp

def calculateCCT():
    from HelpFile import read_input
    data = read_input()

    database = "TCFE12"
    kineticdatabase = "TCFE12"
    dependentmat = data['Material']["Dependentmat"]
    composition = data['Material']["Compositon"]
    phases = ["FCC_A1", "FCC_A1#2", "GAS", "GRAPHITE_A9"]
    dormantphases = ["GAS", "GRAPHITE_A9"]
    referencestates = {"C": "Graphite_A9", "N": "GAS"}

    with TCPython() as start:
        # create and configure a single equilibrium calculation
        calculation = (
            start
            .select_thermodynamic_and_kinetic_databases_with_elements(database, kineticdatabase, [dependentmat] + list(composition))
            .get_system()
            .with_property_model_calculation('TTT Diagram')
            .set_condition()
            .set_phase_to_suspended('*')
            .select_phase("FCC_A1")
            #.disable_global_minimization()
        )
        TTT = calculation.with_ttt_precipitation_calculation()
            #TTT.set_composition()
        ttt_results = (calculation.with_ttt_precipitation_calculation()
                       .set_composition_unit(CompositionUnit.MASS_FRACTION)
                       .set_composition("C", 10)
                       .set_composition("N", 10)
                       .with_matrix_phase(matrix)
                       .set_min_temperature(1000)
                       .set_max_temperature(1160)
                       .set_temperature_step(10)
                       .set_max_annealing_time(1.0e6)
                       .stop_at_volume_fraction_of_phase(1.e-4)
                       .calculate()
                       )

        for element in composition:
            calculation.set_condition(ThermodynamicQuantity.mass_fraction_of_a_component(element), composition[element]/100)
        for phase in phases:
            calculation.set_phase_to_entered(phase)
        for phase in dormantphases:
            calculation.set_phase_to_dormant(phase)
        for element in referencestates:
            calculation.with_reference_state(element, referencestates[element])
        calc_result = (calculation
                       .calculate()  # Aktiverar beräkningen
                       )
        return
def calculatePerlite(temperatures, composition):
    print("Perlite model")
    from HelpFile import read_input
    data = read_input()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = data['Material']["Dependentmat"]
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
    from HelpFile import read_input
    data = read_input()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = data['Material']["Dependentmat"]
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
    from HelpFile import read_input
    data = read_input()
    database = "TCFE12"
    kindatabase = "MOBFE7"
    dependentmat = data['Material']["Dependentmat"]
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

        print("Available arguments: {}".format(calculation.get_arguments()))
        calc_result = (calculation.calculate()  # Aktiverar beräkningen
                       )
        print("Available result quantities: {}".format(calc_result.get_result_quantities()))
        start = calc_result.get_value_of("Ms")
        half = calc_result.get_value_of("M50")
        finish = calc_result.get_value_of("M99")
        print([start, half, finish])
    return start, half, finish
def yieldStrength():
    # Yield Strength
    pass