from simulation import SimRecord, single_tank_simulate, pressurised_liquid_simulate
from cases import n2o_blowdown_case, n2_blowdown_case, ethanol_case

from matplotlib import pyplot as plt


def log_results(file_path:str, sim_record: SimRecord):
    """
    saves the simulation results to a csv file for later analysis, with data rounded to 2 decimal places for readability.
    Works dynamically with any number/name of tanks, injectors, and regulators.
    """
    with open(file_path, 'w') as f:
        # Build header dynamically based on available data
        headers = ["time_s"]
        
        # Add tank headers
        if sim_record.points and sim_record.points[0].tanks is not None:
            for tank_name in sim_record.points[0].tanks.keys():
                headers.extend([
                    f"{tank_name}_pressure_bar",
                    f"{tank_name}_temperature_k",
                    f"{tank_name}_liquid_mass_kg",
                    f"{tank_name}_vapour_mass_kg",
                    f"{tank_name}_pressurant_gas_mass_kg",
                    f"{tank_name}_total_mass_kg",
                    f"{tank_name}_total_internal_energy_j"
                ])
        
        # Add injector headers
        if sim_record.points and sim_record.points[0].injectors_mdot is not None:
            for injector_name in sim_record.points[0].injectors_mdot.keys():
                headers.append(f"{injector_name}_mdot_kg_s")
        
        f.write(",".join(headers) + "\n")
        
        # Write data rows
        for point in sim_record.points:
            row_data = [f"{point.time_s:.2f}"]
            
            # Extract tank data
            if point.tanks is not None:
                for tank_name, tank in point.tanks.items():
                    # Safely format values, writing empty string for None
                    row_data.extend([
                        f"{tank.pressure_pa/1e5:.2f}" if getattr(tank, 'pressure_pa', None) is not None else "",
                        f"{tank.temperature_k:.2f}" if getattr(tank, 'temperature_k', None) is not None else "",
                        f"{tank.liquid_mass_kg:.2f}" if getattr(tank, 'liquid_mass_kg', None) is not None else "",
                        f"{tank.vapour_mass_kg:.2f}" if getattr(tank, 'vapour_mass_kg', None) is not None else "",
                        f"{tank.pressurant_gas_mass_kg:.2f}" if getattr(tank, 'pressurant_gas_mass_kg', None) is not None else "",
                        f"{tank.total_mass_kg:.2f}" if getattr(tank, 'total_mass_kg', None) is not None else "",
                        f"{tank.total_internal_energy_j:.2f}" if getattr(tank, 'total_internal_energy_j', None) is not None else ""
                    ])
            
            # Extract injector data
            if point.injectors_mdot is not None:
                for injector_name, mdot in point.injectors_mdot.items():
                    row_data.append(f"{mdot:.2f}")
            
            f.write(",".join(row_data) + "\n")
    

if __name__ == "__main__":

    """
    run case and plot results
    """
    case = ethanol_case

    sim_record = pressurised_liquid_simulate(case, record=True)

    if sim_record is None:
        raise ValueError("Simulation completed without recording.")
    else:
        print(f"Simulation completed with {len(sim_record.points)} recorded points.")
    
    # build dynamic series dict: keys -> list of values over time
    series = {}
    time_s_list = []

    # determine names from first recorded point (if available)
    first = sim_record.points[0] if sim_record.points else None
    tank_names = list(first.tanks.keys()) if first and first.tanks else []
    injector_names = list(first.injectors_mdot.keys()) if first and first.injectors_mdot else []
    regulator_names = list(first.injectors_mdot.keys()) if first and first.injectors_mdot else []

    # prepare keys for tanks
    for t in tank_names:
        series[f"{t}_pressure_bar"] = []
        series[f"{t}_temperature_k"] = []
        series[f"{t}_liquid_mass_kg"] = []
        series[f"{t}_vapour_mass_kg"] = []
        series[f"{t}_pressurant_gas_mass_kg"] = []
        series[f"{t}_total_mass_kg"] = []
        series[f"{t}_total_internal_energy_kj"] = []

    # prepare keys for injectors
    for inj in injector_names:
        series[f"{inj}_mdot_kg_s"] = []

    # populate series
    for point in sim_record.points:
        time_s_list.append(point.time_s)

        # tanks
        if point.tanks is not None:
            # populate combined mass plots
            for t in tank_names:
                tank = point.tanks.get(t)
                first_tank = first.tanks.get(t) if first and first.tanks else None
                model = getattr(first_tank, 'phase_model', None)
            for t in tank_names:
                tank = point.tanks.get(t)
                if tank is None:
                    series[f"{t}_pressure_bar"].append(None)
                    series[f"{t}_temperature_k"].append(None)
                    series[f"{t}_liquid_mass_kg"].append(None)
                    series[f"{t}_vapour_mass_kg"].append(None)
                    series[f"{t}_pressurant_gas_mass_kg"].append(None)
                    series[f"{t}_total_mass_kg"].append(None)
                    series[f"{t}_total_internal_energy_kj"].append(None)
                else:
                    series[f"{t}_pressure_bar"].append(tank.pressure_pa/1e5 if tank.pressure_pa is not None else None)
                    series[f"{t}_temperature_k"].append(tank.temperature_k)
                    series[f"{t}_liquid_mass_kg"].append(tank.liquid_mass_kg)
                    series[f"{t}_vapour_mass_kg"].append(tank.vapour_mass_kg)
                    series[f"{t}_pressurant_gas_mass_kg"].append(tank.pressurant_gas_mass_kg)
                    series[f"{t}_total_mass_kg"].append(tank.total_mass_kg)
                    series[f"{t}_total_internal_energy_kj"].append(tank.total_internal_energy_j/1e3 if tank.total_internal_energy_j is not None else None)
        else:
            for t in tank_names:
                series[f"{t}_pressure_bar"].append(None)
                series[f"{t}_temperature_k"].append(None)
                series[f"{t}_liquid_mass_kg"].append(None)
                series[f"{t}_vapour_mass_kg"].append(None)
                series[f"{t}_pressurant_gas_mass_kg"].append(None)
                series[f"{t}_total_mass_kg"].append(None)
                series[f"{t}_total_internal_energy_kj"].append(None)

        # injectors
        if point.injectors_mdot is not None:
            for inj in injector_names:
                series[f"{inj}_mdot_kg_s"].append(point.injectors_mdot.get(inj))
        else:
            for inj in injector_names:
                series[f"{inj}_mdot_kg_s"].append(None)

   # create plots dynamically
    keys = list(series.keys())

    # add combined mass plots for tanks that are self_pressurised or pressurised_liquid
    mass_models = {"self_pressurised", "pressurised_liquid"}
    for t in tank_names:
        first_tank = first.tanks.get(t) if first and first.tanks else None
        model = getattr(first_tank.config, 'phase_model', None) if first_tank is not None else None
        if model in mass_models:
            keys.append(f"{t}_all_masses")

    n = len(keys)
    if n == 0:
        print("No series to plot.")
    else:
        cols = 3
        rows = (n + cols - 1) // cols
        fig, axs = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows))
        axs_flat = axs.flatten() if hasattr(axs, 'flatten') else [axs]

        i = -1
        for i, k in enumerate(keys):
            ax = axs_flat[i]

            if k.endswith("_all_masses"):
                t = k[:-11]

                ax.plot(time_s_list, series[f"{t}_liquid_mass_kg"], label="liquid mass")
                ax.plot(time_s_list, series[f"{t}_vapour_mass_kg"], label="vapour mass")
                ax.plot(time_s_list, series[f"{t}_pressurant_gas_mass_kg"], label="pressurant mass")
                ax.plot(time_s_list, series[f"{t}_total_mass_kg"], label="total mass")

                ax.set_title(f"{t} all masses")
                ax.set_xlabel("Time (s)")
                ax.set_ylabel("Mass (kg)")
                ax.legend()

            else:
                ax.plot(time_s_list, series[k])
                ax.set_title(k.replace("_", " "))
                ax.set_xlabel("Time (s)")

                if "pressure" in k:
                    ax.set_ylabel("Pressure (bar)")
                    ax.set_ylim(bottom=0)
                elif "temperature" in k:
                    ax.set_ylabel("Temperature (K)")
                elif "mass" in k:
                    ax.set_ylabel("Mass (kg)")
                    ax.set_ylim(bottom=0)
                elif "energy" in k:
                    ax.set_ylabel("Energy (kJ)")
                elif "mdot" in k:
                    ax.set_ylabel("Mass flow rate (kg/s)")
                    ax.set_ylim(bottom=0)

        # hide unused subplots
        for j in range(i + 1, len(axs_flat)):
            axs_flat[j].axis("off")

        plt.tight_layout()
        fig.savefig(f"{case.name}_test_results.png", dpi=300)
        plt.show()

        log_results(f"{case.name}_results.csv", sim_record)