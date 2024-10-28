from tc_python import *

def TC_Dictra_Err():
    nrLoops = 20
    composition = {"C": 0.2, "Cr": 1.6, "Mn": 0.5, "Ni": 1.5, "Mo": 0.3}
    boost_t = 1
    diff_t = 1
    total_t = (boost_t + diff_t) * nrLoops

    with TCPython() as session:
        logging.getLogger("tc_python").setLevel(logging.INFO)

        TCsystem = (session
                  .select_thermodynamic_and_kinetic_databases_with_elements("TCFE12", "MOBFE7",
                                                                            ["Fe"] + list(composition.keys()))
                  .without_default_phases()
                  .select_phase("FCC_A1")
                  .select_phase("FCC_A1#2")
                  .select_phase("CEMENTITE_D011")
                  .select_phase("GRAPHITE_A9")
                  .get_system())

        compprofile = CompositionProfile(Unit.MASS_PERCENT)
        for element in composition:
            compprofile.add(element, ElementProfile.constant(composition[element]))

        austenite = (Region("Austenite")
                     .set_width(0.008)
                     .with_grid(CalculatedGrid.geometric()
                            .set_no_of_points(50)
                            .set_geometrical_factor(0.9))
                     .add_phase("FCC_A1")
                     .add_phase("FCC_A1#2")
                     .add_phase("CEMENTITE_D011")
                     .with_composition_profile(compprofile))

        calculation = (TCsystem
                       .with_isothermal_diffusion_calculation()
                       .with_solver(Solver.homogenization().with_function(HomogenizationFunctions.labyrinth_factor_f2('FCC_A1#1')))
                       .with_reference_state("C", "GRAPHITE_A9")
                       .set_temperature(1173.15)
                       .set_simulation_time(total_t)
                       .with_cylindrical_geometry()
                       .remove_all_regions()
                       .add_region(austenite)
                       .with_right_boundary_condition(BoundaryCondition.mixed_zero_flux_and_activity()
                                                      .set_activity_for_element("C", str(1.0)))
                       .set_simulation_time(boost_t))

        calculation = calculation.calculate()

        # First diffusion
        calculation = (calculation.with_continued_calculation()
                       .with_right_boundary_condition(BoundaryCondition.closed_system())
                       .set_simulation_time(boost_t + diff_t))
        calculation = calculation.calculate()

        # Looping boost and diffusions
        for i in range(nrLoops):
            # Boost cycle
            calculation = (calculation.with_continued_calculation()
                           .with_right_boundary_condition(BoundaryCondition.mixed_zero_flux_and_activity()
                                                          .set_activity_for_element("C", str(1.0)))
                           .set_simulation_time(boost_t*(i+1)+diff_t*i))
            calculation = calculation.calculate()
            # Diffusion cycle
            calculation = (calculation.with_continued_calculation()
                           .with_right_boundary_condition(BoundaryCondition.closed_system())
                           .set_simulation_time(boost_t*(i+1) + diff_t*(i+1)))
            calculation = calculation.calculate()

        # Extracting results
        distance, mass_frac_C = calculation.get_mass_fraction_of_component_at_time("C", SimulationTime.LAST)
        print(mass_frac_C)
