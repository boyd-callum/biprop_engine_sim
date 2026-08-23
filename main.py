from simulation import biprop_simulate
from cases import full_biprop_case, squirtle_case

from outputs import plot_sim_record, log_results, print_sim_summary
from sweep import *

import numpy as np
    

if __name__ == "__main__":

    """
    run case and plot results
    """
    # case = squirtle_case

    # sim_record = biprop_simulate(case, record=True)

    # if sim_record is None:
    #     raise ValueError("Simulation completed without recording.")
    # else:
    #     print(f"Simulation completed with {len(sim_record.points)} recorded points.")
    
    # print_sim_summary(
    #     simRecord=sim_record,
    #     liquidTankName="n2o_tank"
    # )

    # plot_sim_record(
    #     simRecord=sim_record,
    #     file_path=f"plots/{case.name}_test_results.png",
    #     show=False
    # )

    # log_results(
    #     file_path=f"logs/{case.name}_results.csv",
    #     sim_record=sim_record
    # )


    # --------------------------------------------------------
    # N2O injector area
    # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=range(12, 26),
    #     modify_case=set_n2o_injector_area,
    #     parameter_name="n2o_injector_area_mm2"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="n2o_injector_area_mm2",
    #     file_prefix="squirtle_n2o_injector_area"
    # )


    # # --------------------------------------------------------
    # # Ethanol injector area
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(2.0, 4.1, 0.2),
    #     modify_case=set_ethanol_injector_area,
    #     parameter_name="ethanol_injector_area_mm2"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="ethanol_injector_area_mm2",
    #     file_prefix="squirtle_ethanol_injector_area"
    # )


    # # --------------------------------------------------------
    # # Throat diameter
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(18.0, 29.0, 1.0),
    #     modify_case=set_throat_diameter,
    #     parameter_name="throat_diameter_mm"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="throat_diameter_mm",
    #     file_prefix="squirtle_throat_diameter"
    # )


    # # --------------------------------------------------------
    # # Ethanol regulator pressure
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(35.0, 56.0, 2.5),
    #     modify_case=set_ethanol_pressure,
    #     parameter_name="ethanol_pressure_bar"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="ethanol_pressure_bar",
    #     file_prefix="squirtle_ethanol_pressure"
    # )


    # # --------------------------------------------------------
    # # Initial N2O pressure
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(45.0, 66.0, 2.5),
    #     modify_case=set_n2o_initial_pressure,
    #     parameter_name="n2o_initial_pressure_bar"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="n2o_initial_pressure_bar",
    #     file_prefix="squirtle_n2o_initial_pressure"
    # )


    # # --------------------------------------------------------
    # # N2O mass
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(3.0, 6.1, 0.25),
    #     modify_case=set_n2o_mass,
    #     parameter_name="n2o_mass_kg"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="n2o_mass_kg",
    #     file_prefix="squirtle_n2o_mass"
    # )


    # # --------------------------------------------------------
    # # Ethanol mass
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(1.0, 2.01, 0.1),
    #     modify_case=set_ethanol_mass,
    #     parameter_name="ethanol_mass_kg"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="ethanol_mass_kg",
    #     file_prefix="squirtle_ethanol_mass"
    # )



    # # --------------------------------------------------------
    # # Expansion ratio
    # # --------------------------------------------------------

    # results = run_1D_sweep(
    #     squirtle_case,
    #     values=np.arange(2.0, 8.1, 0.5),
    #     modify_case=set_expansion_ratio,
    #     parameter_name="expansion_ratio"
    # )

    # plot_standard_1D_sweep_outputs(
    #     results,
    #     x_parameter="expansion_ratio",
    #     file_prefix="squirtle_expansion_ratio"
    # )

    # # --------------------------------------------------------
    # # Timestep convergence
    # # --------------------------------------------------------

    results = run_1D_sweep(
        squirtle_case,
        values=[
            0.100,
            0.050,
            0.025,
            0.020,
            0.010,
            0.005,
            0.0025,
            0.001,
        ],
        modify_case=set_timestep,
        parameter_name="dt_s"
    )

    plot_1D_sweep(
        results,
        x_parameter="dt_s",
        y_parameter="total_impulse_ns",
        file_path="plots/squirtle_timestep_vs_impulse.png"
    )

    plot_1D_sweep(
        results,
        x_parameter="dt_s",
        y_parameter="liquid_avg_of",
        file_path="plots/squirtle_timestep_vs_liquid_of.png"
    )

    plot_1D_sweep(
        results,
        x_parameter="dt_s",
        y_parameter="liquid_avg_thrust_n",
        file_path="plots/squirtle_timestep_vs_liquid_thrust.png"
    )

    plot_1D_sweep(
        results,
        x_parameter="dt_s",
        y_parameter="liquid_avg_chamber_pressure_bar",
        file_path="plots/squirtle_timestep_vs_liquid_pc.png"
    )