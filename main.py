from simulation import SimRecord, blowdown_simulate
from cases import n2o_blowdown_case

from matplotlib import pyplot as plt


def log_results(file_path:str, sim_record: SimRecord):
    """
    saves the simulation results to a csv file for later analysis, with data rounded to 2 decimal places for readability.
    """
    with open(file_path, 'w') as f:
        f.write("time_s,pressure_bar,temperature_k,liquid_mass_kg,vapour_mass_kg,total_mass_kg,injector_mdot_kg_s,total_internal_energy_j\n")
        for point in sim_record.points:
            time_s = point.time_s
            tank_pressure_pa = point.tanks['n2o_tank'].pressure_pa if point.tanks is not None else None
            tank_temperature_k = point.tanks['n2o_tank'].temperature_k if point.tanks is not None else None
            liquid_mass_kg = point.tanks['n2o_tank'].liquid_mass_kg if point.tanks is not None else None
            vapour_mass_kg = point.tanks['n2o_tank'].vapour_mass_kg if point.tanks is not None else None
            total_mass_kg = point.tanks['n2o_tank'].total_mass_kg if point.tanks is not None else None
            mdot_kg_s = point.injectors_mdot['test_injector'] if point.injectors_mdot is not None else None
            total_internal_energy_j = point.tanks['n2o_tank'].total_internal_energy_j if point.tanks is not None else None

            if tank_pressure_pa is not None:
                f.write(f"{time_s:.2f},{tank_pressure_pa/1e5:.2f},{tank_temperature_k:.2f},{liquid_mass_kg:.2f},{vapour_mass_kg:.2f},{total_mass_kg:.2f},{mdot_kg_s:.2f},{total_internal_energy_j:.2f}\n")
    

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
    total_energy_j_list = []

    for point in sim_record.points:
        time_s = point.time_s
        tank_pressure_pa = point.tanks['n2o_tank'].pressure_pa if point.tanks is not None else None
        tank_temperature_k = point.tanks['n2o_tank'].temperature_k if point.tanks is not None else None
        liquid_mass_kg = point.tanks['n2o_tank'].liquid_mass_kg if point.tanks is not None else None
        vapour_mass_kg = point.tanks['n2o_tank'].vapour_mass_kg if point.tanks is not None else None
        total_mass_kg = point.tanks['n2o_tank'].total_mass_kg if point.tanks is not None else None
        mdot_kg_s = point.injectors_mdot['test_injector'] if point.injectors_mdot is not None else None
        total_energy_j = point.tanks['n2o_tank'].total_internal_energy_j if point.tanks is not None else None
    
        time_s_list.append(time_s)
        pressure_pa_list.append(tank_pressure_pa)
        temperature_k_list.append(tank_temperature_k)
        liquid_mass_kg_list.append(liquid_mass_kg)
        vapour_mass_kg_list.append(vapour_mass_kg)
        total_mass_kg_list.append(total_mass_kg)
        mdot_kg_s_list.append(mdot_kg_s)
        total_energy_j_list.append(total_energy_j)
    
    # plot results
    # pressure, temp, liauid mass, vapour mass, total mass, mdot, (liquid, vapour and total) mass on same plot, total energy
    fig, axs = plt.subplots(3, 3, figsize=(12, 10))
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
    axs[0, 2].plot(time_s_list, liquid_mass_kg_list, label="Liquid Mass")
    axs[0, 2].plot(time_s_list, vapour_mass_kg_list, label="Vapour Mass")
    axs[0, 2].plot(time_s_list, total_mass_kg_list, label="Total Mass")
    axs[0, 2].set_title("Masses (kg)")
    axs[0, 2].legend()
    axs[1, 2].plot(time_s_list, total_energy_j_list)
    axs[1, 2].set_title("Total Internal Energy (J)")
    plt.tight_layout()
    plt.show()


    log_results("blowdown_test_results.csv", sim_record)