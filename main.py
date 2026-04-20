from simulation import blowdown_simulate
from cases import n2o_blowdown_case

from matplotlib import pyplot as plt


if __name__ == "__main__":

    """
    run case and plot results
    """

    sim_record = blowdown_simulate(n2o_blowdown_case, record=True)

    if sim_record is None:
        raise ValueError("Simulation completed without recording.")
    else:
        print(f"Simulation completed with {len(sim_record.points)} recorded points.")
    
    # set up lists to hold time and tank data for plotting
    time_s_list = []
    pressure_pa_list = []
    temperature_k_list = []
    liquid_mass_kg_list = []
    vapour_mass_kg_list = []
    total_mass_kg_list = []
    mdot_kg_s_list = []

    for point in sim_record.points:
        time_s = point.time_s
        tank_pressure_pa = point.tanks['n2o_tank'].pressure_pa if point.tanks is not None else None
        tank_temperature_k = point.tanks['n2o_tank'].temperature_k if point.tanks is not None else None
        liquid_mass_kg = point.tanks['n2o_tank'].liquid_mass_kg if point.tanks is not None else None
        vapour_mass_kg = point.tanks['n2o_tank'].vapour_mass_kg if point.tanks is not None else None
        total_mass_kg = point.tanks['n2o_tank'].total_mass_kg if point.tanks is not None else None
        mdot_kg_s = point.injectors_mdot['test_injector'] if point.injectors_mdot is not None else None
    
        time_s_list.append(time_s)
        pressure_pa_list.append(tank_pressure_pa)
        temperature_k_list.append(tank_temperature_k)
        liquid_mass_kg_list.append(liquid_mass_kg)
        vapour_mass_kg_list.append(vapour_mass_kg)
        total_mass_kg_list.append(total_mass_kg)
        mdot_kg_s_list.append(mdot_kg_s)
    
    # plot results
    fig, axs = plt.subplots(3, 2, figsize=(12, 10))
    axs[0, 0].plot(time_s_list, pressure_pa_list)
    axs[0, 0].set_title("Tank Pressure (Pa)")
    axs[0, 1].plot(time_s_list, temperature_k_list)
    axs[0, 1].set_title("Tank Temperature (K)")
    axs[1, 0].plot(time_s_list, liquid_mass_kg_list)
    axs[1, 0].set_title("Liquid Mass (kg)")
    axs[1, 1].plot(time_s_list, vapour_mass_kg_list)
    axs[1, 1].set_title("Vapour Mass (kg)")
    axs[2, 0].plot(time_s_list, total_mass_kg_list)
    axs[2, 0].set_title("Total Mass (kg)")
    axs[2, 1].plot(time_s_list, mdot_kg_s_list)
    axs[2, 1].set_title("Injector Mass Flow Rate (kg/s)")
    plt.tight_layout()
    plt.show()


    