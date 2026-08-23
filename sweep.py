from dataclasses import dataclass, replace
import pandas as pd
from typing import Callable
from pathlib import Path
import matplotlib.pyplot as plt
import math

from cases import SimCase

from tanks import TankConfig, TankInitialCondition
from fluid import Fluid
from injector import InjectorConfig
from engine import EngineConfig, EngineGeometry
from regulator import RegulatorConfig

from simulation import biprop_simulate
from outputs import calculate_sim_summary, SimSummary




def run_1D_sweep(
    base_case: SimCase,
    values: list[float],
    modify_case: Callable,
    parameter_name: str
) -> pd.dataframe:

    results = []

    for value in values:

        case = modify_case(base_case, value)

        sim_record = biprop_simulate(case)

        summary = calculate_sim_summary(sim_record, liquidTankName="n2o_tank")

        results.append({
            parameter_name: value,

            "avg_of": summary.avg_of,
            "liquid_avg_of": summary.liquid_avg_of,

            "avg_thrust_n": summary.avg_thrust_n,
            "liquid_avg_thrust_n": summary.liquid_avg_thrust_n,
            "peak_thrust_n": summary.peak_thrust_n,

            "avg_isp_s": summary.avg_isp_s,
            "liquid_avg_isp_s": summary.liquid_avg_isp_s,

            "total_impulse_ns": summary.total_impulse_ns,

            "avg_chamber_pressure_bar":
                summary.avg_chamber_pressure_bar,
            "liquid_avg_chamber_pressure_bar": summary.liquid_avg_chamber_pressure_bar,

            "peak_chamber_pressure_bar":
                summary.peak_chamber_pressure_bar,
        })

        print(
            f"{parameter_name}: {value} | "
            f"Impulse: {summary.total_impulse_ns:.1f} Ns | "
            f"Liquid O/F: {summary.liquid_avg_of:.3f} | "
            f"Liquid thrust: {summary.liquid_avg_thrust_n:.1f} N"
        )

    return pd.DataFrame(results)


def plot_1D_sweep(
        results: pd.DataFrame,
        x_parameter: str,
        y_parameter: str,
        file_path: str | Path | None = None,
        show: bool = False
) -> None:

    fig = plt.figure(figsize=(8,5))

    plt.plot(
        results[x_parameter],
        results[y_parameter]
    )
    plt.xlabel(x_parameter)
    plt.ylabel(y_parameter)
    plt.title(f"{y_parameter} vs {x_parameter}")

    plt.grid(True)
    plt.tight_layout()
    
    if file_path is not None:
        fig.savefig(
            file_path,
            dpi=300,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    plt.close(fig)

def plot_standard_1D_sweep_outputs(
    results,
    x_parameter: str,
    file_prefix: str,
) -> None:

    plot_1D_sweep(
        results,
        x_parameter=x_parameter,
        y_parameter="total_impulse_ns",
        file_path=f"plots/{file_prefix}_vs_impulse.png"
    )

    plot_1D_sweep(
        results,
        x_parameter=x_parameter,
        y_parameter="liquid_avg_of",
        file_path=f"plots/{file_prefix}_vs_liquid_of.png"
    )

    plot_1D_sweep(
        results,
        x_parameter=x_parameter,
        y_parameter="liquid_avg_thrust_n",
        file_path=f"plots/{file_prefix}_vs_liquid_thrust.png"
    )

    plot_1D_sweep(
        results,
        x_parameter=x_parameter,
        y_parameter="liquid_avg_chamber_pressure_bar",
        file_path=f"plots/{file_prefix}_vs_liquid_pc.png"
    )

def set_n2o_injector_area(
    case: SimCase,
    area_mm2: float
) -> SimCase:



    injector = replace(
        case.injector_configs["n2o_injector"],
        area_m2 = area_mm2*1e-6 # convert from mm2 to m2
        )

    new_case = replace(
        case,
        injector_configs = {
            **case.injector_configs,
            "n2o_injector": injector
        }
    )

    return new_case


def set_ethanol_injector_area(
    case: SimCase,
    area_mm2: float,
) -> SimCase:

    injector = replace(
        case.injector_configs["ethanol_injector"],
        area_m2=area_mm2 * 1e-6,
    )

    return replace(
        case,
        injector_configs={
            **case.injector_configs,
            "ethanol_injector": injector,
        },
    )

def set_throat_diameter(
    case: SimCase,
    diameter_mm: float,
) -> SimCase:

    diameter_m = diameter_mm * 1e-3
    throat_area_m2 = math.pi * diameter_m**2 / 4

    geometry = replace(
        case.engine_config.geometry,
        nozzle_throat_area_m2=throat_area_m2,
    )

    engine = replace(
        case.engine_config,
        geometry=geometry,
    )

    return replace(
        case,
        engine_config=engine,
    )

def set_ethanol_pressure(
    case: SimCase,
    pressure_bar: float,
) -> SimCase:

    pressure_pa = pressure_bar * 1e5

    regulator = replace(
        case.regulator_configs["ethanol_regulator"],
        set_pressure_pa=pressure_pa,
    )

    initial = replace(
        case.tank_initial_conditions["ethanol_tank"],
        pressure_pa=pressure_pa,
    )

    return replace(
        case,

        regulator_configs={
            **case.regulator_configs,
            "ethanol_regulator": regulator,
        },

        tank_initial_conditions={
            **case.tank_initial_conditions,
            "ethanol_tank": initial,
        },
    )

def set_n2o_initial_pressure(
    case: SimCase,
    pressure_bar: float,
) -> SimCase:

    initial = replace(
        case.tank_initial_conditions["n2o_tank"],
        pressure_pa=pressure_bar * 1e5,
    )

    return replace(
        case,
        tank_initial_conditions={
            **case.tank_initial_conditions,
            "n2o_tank": initial,
        },
    )

def set_n2o_mass(
    case: SimCase,
    mass_kg: float,
) -> SimCase:

    initial = replace(
        case.tank_initial_conditions["n2o_tank"],
        total_mass_kg=mass_kg,
    )

    return replace(
        case,
        tank_initial_conditions={
            **case.tank_initial_conditions,
            "n2o_tank": initial,
        },
    )

def set_ethanol_mass(
    case: SimCase,
    mass_kg: float,
) -> SimCase:

    initial = replace(
        case.tank_initial_conditions["ethanol_tank"],
        total_mass_kg=mass_kg,
    )

    return replace(
        case,
        tank_initial_conditions={
            **case.tank_initial_conditions,
            "ethanol_tank": initial,
        },
    )

def set_expansion_ratio(
    case: SimCase,
    expansion_ratio: float,
) -> SimCase:

    geometry = replace(
        case.engine_config.geometry,
        expansion_ratio=expansion_ratio,
    )

    engine = replace(
        case.engine_config,
        geometry=geometry,
    )

    return replace(
        case,
        engine_config=engine,
    )

def set_timestep(
    case: SimCase,
    dt_s: float,
) -> SimCase:

    settings = replace(
        case.settings,
        dt_s=dt_s,
    )

    return replace(
        case,
        settings=settings,
    )