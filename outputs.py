from __future__ import annotations
from pathlib import Path
from typing import Any, Iterable, TYPE_CHECKING
import math
from dataclasses import dataclass

import matplotlib.pyplot as plt

import csv

from simulation import SimRecord, SimPoint
from helpers import safe_divide, get_single_config_by_role, get_tank_phase_model

if TYPE_CHECKING:
    from tanks import PhaseModel


def append_tank_series(
    series: dict[str, list[float | None]],
    tankName: str,
    tankState: Any | None,
) -> None:
    """
    Appends one timestep of tank data to the series dictionary.
    """

    if tankState is None:
        series[f"{tankName}_pressure_bar"].append(None)
        series[f"{tankName}_temperature_k"].append(None)
        series[f"{tankName}_liquid_mass_kg"].append(None)
        series[f"{tankName}_vapour_mass_kg"].append(None)
        series[f"{tankName}_pressurant_gas_mass_kg"].append(None)
        series[f"{tankName}_total_mass_kg"].append(None)
        series[f"{tankName}_total_internal_energy_kj"].append(None)
        return

    series[f"{tankName}_pressure_bar"].append(
        safe_divide(tankState.pressure_pa, 1e5)
    )
    series[f"{tankName}_temperature_k"].append(tankState.temperature_k)
    series[f"{tankName}_liquid_mass_kg"].append(tankState.liquid_mass_kg)
    series[f"{tankName}_vapour_mass_kg"].append(tankState.vapour_mass_kg)
    series[f"{tankName}_pressurant_gas_mass_kg"].append(
        tankState.pressurant_gas_mass_kg
    )
    series[f"{tankName}_total_mass_kg"].append(tankState.total_mass_kg)
    series[f"{tankName}_total_internal_energy_kj"].append(
        safe_divide(tankState.total_internal_energy_j, 1e3)
    )


def append_engine_series(
    series: dict[str, list[float | None]],
    engineState: Any | None,
) -> None:
    """
    Appends one timestep of engine data to the series dictionary.
    """

    if engineState is None:
        series["engine_chamber_pressure_bar"].append(None)
        series["engine_chamber_temperature_k"].append(None)
        series["engine_thrust_n"].append(None)
        series["engine_isp_s"].append(None)
        series["engine_mixture_ratio"].append(None)
        series["engine_total_mdot_kg_s"].append(None)
        series["engine_cstar_ideal_m_s"].append(None)
        series["engine_cstar_delivered_m_s"].append(None)
        series["engine_cf_ideal"].append(None)
        series["engine_cf_delivered"].append(None)
        return

    series["engine_chamber_pressure_bar"].append(
        safe_divide(engineState.chamber_pressure_pa, 1e5)
    )
    series["engine_chamber_temperature_k"].append(
        engineState.chamber_temperature_k
    )
    series["engine_thrust_n"].append(engineState.thrust_n)
    series["engine_isp_s"].append(engineState.isp_s)
    series["engine_mixture_ratio"].append(engineState.mixture_ratio)
    series["engine_total_mdot_kg_s"].append(engineState.total_mdot_kg_s)
    series["engine_cstar_ideal_m_s"].append(engineState.cstar_ideal_m_s)
    series["engine_cstar_delivered_m_s"].append(
        engineState.cstar_delivered_m_s
    )
    series["engine_cf_ideal"].append(engineState.cf_ideal)
    series["engine_cf_delivered"].append(engineState.cf_delivered)



def build_sim_series(
    simRecord: SimRecord,
) -> tuple[
    list[float],
    dict[str, list[float | None]],
    list[str],
    list[str],
    list[str],
]:
    """
    Builds plottable time-series data from a SimRecord.

    Returns:
        tuple:
            ``(timeSList, series, tankNames, injectorNames, regulatorNames)``.
    """

    if len(simRecord.points) == 0:
        return [], {}, [], [], []

    firstPoint = simRecord.points[0]

    tankNames = list(firstPoint.tanks.keys()) if firstPoint.tanks else []
    injectorNames = (
        list(firstPoint.injectors_mdot.keys())
        if firstPoint.injectors_mdot
        else []
    )
    regulatorNames = (
        list(firstPoint.regulators_mdot.keys())
        if firstPoint.regulators_mdot
        else []
    )

    timeSList: list[float] = []
    series: dict[str, list[float | None]] = {}

    for tankName in tankNames:
        series[f"{tankName}_pressure_bar"] = []
        series[f"{tankName}_temperature_k"] = []
        series[f"{tankName}_liquid_mass_kg"] = []
        series[f"{tankName}_vapour_mass_kg"] = []
        series[f"{tankName}_pressurant_gas_mass_kg"] = []
        series[f"{tankName}_total_mass_kg"] = []
        series[f"{tankName}_total_internal_energy_kj"] = []

    for injectorName in injectorNames:
        series[f"{injectorName}_mdot_kg_s"] = []

    for regulatorName in regulatorNames:
        series[f"{regulatorName}_mdot_kg_s"] = []

    series["engine_chamber_pressure_bar"] = []
    series["engine_chamber_temperature_k"] = []
    series["engine_thrust_n"] = []
    series["engine_isp_s"] = []
    series["engine_mixture_ratio"] = []
    series["engine_total_mdot_kg_s"] = []
    series["engine_cstar_ideal_m_s"] = []
    series["engine_cstar_delivered_m_s"] = []
    series["engine_cf_ideal"] = []
    series["engine_cf_delivered"] = []

    for point in simRecord.points:
        timeSList.append(point.time_s)

        for tankName in tankNames:
            tankState = point.tanks.get(tankName) if point.tanks else None
            append_tank_series(series, tankName, tankState)

        for injectorName in injectorNames:
            mdot = (
                point.injectors_mdot.get(injectorName)
                if point.injectors_mdot
                else None
            )
            series[f"{injectorName}_mdot_kg_s"].append(mdot)

        for regulatorName in regulatorNames:
            mdot = (
                point.regulators_mdot.get(regulatorName)
                if point.regulators_mdot
                else None
            )
            series[f"{regulatorName}_mdot_kg_s"].append(mdot)

        append_engine_series(series, point.engine)

    return timeSList, series, tankNames, injectorNames, regulatorNames


def create_column_axes(
    fig: Any,
    gridCell: Any,
    rowCount: int
) -> list[Any]:

    """
    creates vertically stacked axes within one figure column
    all axes within the column share the same time axis
    """

    subGrid = gridCell.subgridspec(
        rowCount,
        1,
        hspace=0.08
    )

    axes: list[Any] = []

    for rowIndex in range(rowCount):

        if rowIndex == 0:
            ax = fig.add_subplot(subGrid[rowIndex, 0])

        else:
            ax = fig.add_subplot(subGrid[rowIndex, 0], sharex=axes[0])

        axes.append(ax)


    for ax in axes[:-1]:

        ax.tick_params(axis="x", labelbottom=False)

    axes[-1].set_xlabel("Time(s)")

    return axes



def plot_tank_mass(
    ax: Any,
    timeSList: list[float],
    series: dict[str, list[float | None]],
    tankName: str,
    phaseModel: PhaseModel
) -> None:



    """
    plots mass components for each tank model.

    self-pressurised:
        liquid + vapour + total

    pressurised liquid:
        liquid + pressurant + total

    single phase:
        total only
    """


    if phaseModel == "self_pressurised":
        ax.plot(
            timeSList,
            series[f"{tankName}_liquid_mass_kg"],
            label="liquid",
        )

        ax.plot(
            timeSList,
            series[f"{tankName}_vapour_mass_kg"],
            label="vapour",
        )

        ax.plot(
            timeSList,
            series[f"{tankName}_total_mass_kg"],
            label="total",
        )

        ax.legend(fontsize="small")

    elif phaseModel == "pressurised_liquid":
        ax.plot(
            timeSList,
            series[f"{tankName}_liquid_mass_kg"],
            label="liquid",
        )

        ax.plot(
            timeSList,
            series[f"{tankName}_pressurant_gas_mass_kg"],
            label="pressurant",
        )

        ax.plot(
            timeSList,
            series[f"{tankName}_total_mass_kg"],
            label="total",
        )

        ax.legend(fontsize="small")

    else:
        # single_phase
        ax.plot(
            timeSList,
            series[f"{tankName}_total_mass_kg"],
        )

    ax.set_ylim(bottom=0)





def plot_sim_record(
    simRecord: SimRecord,
    file_path: str | Path | None = None,
    show: bool = True,
) -> None:
    """
    Plots simulation results in a single figure.

    Columns:
        - one column per tank
        - feed system
        - engine

    Tank columns:
        - pressure
        - temperature
        - relevant masses
        - internal energy

    Feed-system column:
        - injector mass flows
        - regulator mass flows
        - total engine mass flow

    Engine column:
        - chamber + fuel + oxidiser tank pressures
        - chamber temperature
        - thrust
        - specific impulse
        - mixture ratio
        - characteristic velocity
        - thrust coefficient
    """

    (
        timeSList,
        series,
        tankNames,
        injectorNames,
        regulatorNames,
    ) = build_sim_series(simRecord)

    if len(series) == 0:
        print("No series to plot.")
        return

    firstPoint = simRecord.points[0]

    # --------------------------------------------------------
    # Identify fuel and oxidiser tanks
    # --------------------------------------------------------

    tankConfigs = {
        tankName: tankState.config
        for tankName, tankState in firstPoint.tanks.items()
    }

    oxidiserTankName, _ = get_single_config_by_role(
        tankConfigs,
        "oxidiser",
        "tank",
    )

    fuelTankName, _ = get_single_config_by_role(
        tankConfigs,
        "fuel",
        "tank",
    )

    # --------------------------------------------------------
    # Feed-system plots
    # --------------------------------------------------------

    feedPlots: list[
        tuple[str, str, str]
    ] = []

    for injectorName in injectorNames:
        feedPlots.append(
            (
                injectorName.replace("_", " "),
                f"{injectorName}_mdot_kg_s",
                "Mass flow (kg/s)",
            )
        )

    for regulatorName in regulatorNames:
        feedPlots.append(
            (
                regulatorName.replace("_", " "),
                f"{regulatorName}_mdot_kg_s",
                "Mass flow (kg/s)",
            )
        )

    feedPlots.append(
        (
            "Engine total mass flow",
            "engine_total_mdot_kg_s",
            "Mass flow (kg/s)",
        )
    )

    # --------------------------------------------------------
    # Engine plots
    # --------------------------------------------------------

    enginePlots: list[
        tuple[
            str,
            str,
            list[
                tuple[
                    str,
                    str | None,
                ]
            ],
        ]
    ] = [
        (
            "Pressure",
            "Pressure (bar)",
            [
                (
                    "engine_chamber_pressure_bar",
                    "chamber",
                ),
                (
                    f"{oxidiserTankName}_pressure_bar",
                    "oxidiser tank",
                ),
                (
                    f"{fuelTankName}_pressure_bar",
                    "fuel tank",
                ),
            ],
        ),
        (
            "Chamber temperature",
            "Temperature (K)",
            [
                (
                    "engine_chamber_temperature_k",
                    None,
                ),
            ],
        ),
        (
            "Thrust",
            "Thrust (N)",
            [
                (
                    "engine_thrust_n",
                    None,
                ),
            ],
        ),
        (
            "Specific impulse",
            "Isp (s)",
            [
                (
                    "engine_isp_s",
                    None,
                ),
            ],
        ),
        (
            "Mixture ratio",
            "O/F",
            [
                (
                    "engine_mixture_ratio",
                    None,
                ),
            ],
        ),
        (
            "Characteristic velocity",
            "c* (m/s)",
            [
                (
                    "engine_cstar_ideal_m_s",
                    "ideal",
                ),
                (
                    "engine_cstar_delivered_m_s",
                    "delivered",
                ),
            ],
        ),
        (
            "Thrust coefficient",
            "Cf",
            [
                (
                    "engine_cf_ideal",
                    "ideal",
                ),
                (
                    "engine_cf_delivered",
                    "delivered",
                ),
            ],
        ),
    ]

        # --------------------------------------------------------
    # Main figure
    # --------------------------------------------------------

    totalColumns = len(tankNames) + 2
    maxRows = max(4, len(feedPlots), len(enginePlots))

    fig = plt.figure(
        figsize=(4 * totalColumns, 2.2 * maxRows),
        constrained_layout=True,
    )

    mainGrid = fig.add_gridspec(1, totalColumns, wspace=0.25)

    # --------------------------------------------------------
    # Tank columns
    # --------------------------------------------------------

    for columnIndex, tankName in enumerate(tankNames):
        axes = create_column_axes(fig, mainGrid[0, columnIndex], 4)

        pressureAx = axes[0]
        temperatureAx = axes[1]
        massAx = axes[2]
        energyAx = axes[3]

        tankState = firstPoint.tanks.get(tankName)
        phaseModel = get_tank_phase_model(tankState)
        tankTitle = tankName.replace("_", " ")

        # Pressure
        pressureAx.plot(timeSList, series[f"{tankName}_pressure_bar"])
        pressureAx.set_title(f"{tankTitle} - Pressure")
        pressureAx.set_ylabel("Pressure (bar)")
        pressureAx.set_ylim(bottom=0)

        # Temperature
        temperatureAx.plot(timeSList, series[f"{tankName}_temperature_k"])
        temperatureAx.set_title(f"{tankTitle} - Temperature")
        temperatureAx.set_ylabel("Temperature (K)")

        # Mass
        plot_tank_mass(massAx, timeSList, series, tankName, phaseModel)
        massAx.set_title(f"{tankTitle} - Mass")
        massAx.set_ylabel("Mass (kg)")

        # Internal energy
        energyAx.plot(timeSList, series[f"{tankName}_total_internal_energy_kj"])
        energyAx.set_title(f"{tankTitle} - Internal energy")
        energyAx.set_ylabel("Energy (kJ)")

    # --------------------------------------------------------
    # Feed-system column
    # --------------------------------------------------------

    feedColumnIndex = len(tankNames)

    feedAxes = create_column_axes(
        fig,
        mainGrid[0, feedColumnIndex],
        len(feedPlots),
    )

    for plotIndex, (title, seriesKey, ylabel) in enumerate(feedPlots):
        ax = feedAxes[plotIndex]

        ax.plot(timeSList, series[seriesKey])
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)

    # --------------------------------------------------------
    # Engine column
    # --------------------------------------------------------

    engineColumnIndex = len(tankNames) + 1

    engineAxes = create_column_axes(
        fig,
        mainGrid[0, engineColumnIndex],
        len(enginePlots),
    )

    for plotIndex, (title, ylabel, plotSeries) in enumerate(enginePlots):
        ax = engineAxes[plotIndex]

        for seriesKey, label in plotSeries:
            ax.plot(timeSList, series[seriesKey], label=label)

        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)

        if any(label is not None for _, label in plotSeries):
            ax.legend(fontsize="small")

    # --------------------------------------------------------
    # Save / display
    # --------------------------------------------------------

    if file_path is not None:
        fig.savefig(
            file_path,
            dpi=300,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    plt.close(fig)





def format_csv_value(
    value: float | int | str | None,
    precision: int,
) -> str:
    """
    Formats a value for CSV output, using an empty string for None.
    """

    if value is None:
        return ""

    if isinstance(value, float):
        return f"{value:.{precision}f}"

    return str(value)


def get_ordered_dict_keys(
    dictionaries: Iterable[dict[str, Any] | None],
) -> list[str]:
    """
    Gets unique dictionary keys in the order they first appear.
    """

    orderedKeys: list[str] = []
    seenKeys: set[str] = set()

    for dictionary in dictionaries:
        if dictionary is None:
            continue

        for key in dictionary.keys():
            if key not in seenKeys:
                orderedKeys.append(key)
                seenKeys.add(key)

    return orderedKeys


def log_results(
    file_path: str | Path,
    sim_record: SimRecord,
    precision: int = 6,
) -> None:
    """
    Saves simulation results to a CSV file.

    The output columns are generated dynamically from the recorded simulation
    data, so the function works with any number of tanks, injectors, and
    regulators. Engine outputs are included if engine states are present in
    the simulation record.

    Args:
        file_path:
            Path to the output CSV file.
        sim_record:
            Completed simulation record.
        precision:
            Number of decimal places used for floating-point output.
    """

    filePath = Path(file_path)

    tankNames = get_ordered_dict_keys(
        point.tanks for point in sim_record.points
    )
    injectorNames = get_ordered_dict_keys(
        point.injectors_mdot for point in sim_record.points
    )
    regulatorNames = get_ordered_dict_keys(
        point.regulators_mdot for point in sim_record.points
    )

    enginePresent = any(point.engine is not None for point in sim_record.points)

    tankFields = [
        ("pressure_bar", lambda tank: tank.pressure_pa / 1e5 if tank.pressure_pa is not None else None),
        ("temperature_k", lambda tank: tank.temperature_k),
        ("liquid_mass_kg", lambda tank: tank.liquid_mass_kg),
        ("vapour_mass_kg", lambda tank: tank.vapour_mass_kg),
        ("pressurant_gas_mass_kg", lambda tank: tank.pressurant_gas_mass_kg),
        ("total_mass_kg", lambda tank: tank.total_mass_kg),
        ("total_internal_energy_kj", lambda tank: tank.total_internal_energy_j / 1e3 if tank.total_internal_energy_j is not None else None),
    ]

    engineFields = [
        ("engine_chamber_pressure_bar", lambda engine: engine.chamber_pressure_pa / 1e5 if engine.chamber_pressure_pa is not None else None),
        ("engine_chamber_temperature_k", lambda engine: engine.chamber_temperature_k),
        ("engine_thrust_n", lambda engine: engine.thrust_n),
        ("engine_isp_s", lambda engine: engine.isp_s),
        ("engine_mixture_ratio", lambda engine: engine.mixture_ratio),
        ("engine_total_mdot_kg_s", lambda engine: engine.total_mdot_kg_s),
        ("engine_cstar_ideal_m_s", lambda engine: engine.cstar_ideal_m_s),
        ("engine_cstar_delivered_m_s", lambda engine: engine.cstar_delivered_m_s),
        ("engine_cf_ideal", lambda engine: engine.cf_ideal),
        ("engine_cf_delivered", lambda engine: engine.cf_delivered),
    ]

    headers = ["time_s"]

    for tankName in tankNames:
        for fieldName, _ in tankFields:
            headers.append(f"{tankName}_{fieldName}")

    for injectorName in injectorNames:
        headers.append(f"{injectorName}_mdot_kg_s")

    for regulatorName in regulatorNames:
        headers.append(f"{regulatorName}_mdot_kg_s")

    if enginePresent:
        for fieldName, _ in engineFields:
            headers.append(fieldName)

    with filePath.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headers)

        for point in sim_record.points:
            rowData: list[str] = [
                format_csv_value(point.time_s, precision)
            ]

            for tankName in tankNames:
                tank = point.tanks.get(tankName) if point.tanks is not None else None

                for _, valueGetter in tankFields:
                    value = valueGetter(tank) if tank is not None else None
                    rowData.append(format_csv_value(value, precision))

            for injectorName in injectorNames:
                mdot = (
                    point.injectors_mdot.get(injectorName)
                    if point.injectors_mdot is not None
                    else None
                )
                rowData.append(format_csv_value(mdot, precision))

            for regulatorName in regulatorNames:
                mdot = (
                    point.regulators_mdot.get(regulatorName)
                    if point.regulators_mdot is not None
                    else None
                )
                rowData.append(format_csv_value(mdot, precision))

            if enginePresent:
                engine = point.engine

                for _, valueGetter in engineFields:
                    value = valueGetter(engine) if engine is not None else None
                    rowData.append(format_csv_value(value, precision))

            writer.writerow(rowData)



def get_liquid_phase_mask(
    simRecord: SimRecord,
    tankName: str,
    liquidPhaseThreshold: float = 0.01
) -> list[bool]:

    """
    returns a mask indicating whether the specified tank still has any liquid present.
    """

    liquidMask: list[bool] = []

    for point in simRecord.points:
        tank = point.tanks.get(tankName) if point.tanks is not None else None

        liquidMask.append(tank is not None and tank.liquid_mass_kg is not None and tank.liquid_mass_kg > liquidPhaseThreshold)

    return liquidMask



def get_engine_value(point: SimPoint, attributeName: str) -> float | None:


    if point.engine is None:
        return None

    value = getattr(point.engine, attributeName, None)

    if value is None or not math.isfinite(value):
        return None

    return value

def trap_integrate(
    simRecord: SimRecord,
    attributeName: str,
    mask: list[bool] | None = None
) -> tuple[float, float]:

    integral = 0.0
    duration = 0.0

    for index in range(len(simRecord.points) - 1):
        point0 = simRecord.points[index]
        point1 = simRecord.points[index + 1]

        if mask is not None and not (mask[index] and mask[index + 1]):
            continue

        value0 = get_engine_value(point0, attributeName)
        value1 = get_engine_value(point1, attributeName)

        if value1 is None:
            # nothing known at the end of this interval - the burn has stopped
            continue

        if value0 is None:
            # The recorded initial condition carries no engine state, because the
            # engine is only solved after the tanks have been advanced. The value
            # at the end of the interval was computed from the tank state at its
            # start, so it applies across the whole interval - hold it rather than
            # skipping, which used to drop the first timestep from every integral.
            value0 = value1

        dt = point1.time_s - point0.time_s

        if dt <= 0:
            continue


        integral += 0.5 * (value0+value1) * dt
        duration += dt

    return(integral, duration)


    

def time_average(
    simRecord: SimRecord,
    attributeName: str,
    mask: list[bool] | None = None
) -> float | None:

    integral, duration = trap_integrate(simRecord, attributeName, mask)

    if duration <= 0:
        return None

    return integral / duration


def format_value(
    value: float | None,
    precision: int = 2
) -> str:

    if value is None: 
        return "N/A"

    return f"{value:.{precision}f}"



@dataclass
class SimSummary:
    burn_time_s: float
    liquid_burn_time_s: float

    avg_of: float | None
    liquid_avg_of: float | None

    avg_thrust_n: float | None
    liquid_avg_thrust_n: float | None
    peak_thrust_n: float | None

    avg_chamber_pressure_bar: float | None
    liquid_avg_chamber_pressure_bar: float | None
    peak_chamber_pressure_bar: float | None

    avg_isp_s: float | None
    liquid_avg_isp_s: float | None

    total_impulse_ns: float
    liquid_impulse_ns: float




def calculate_sim_summary(
    simRecord: SimRecord,
    liquidTankName: str,
) -> SimSummary:
    """
    calculates overall and liquid-phase engine performance stats.

    - avg O/F
    - avg thrust
    - peak thrust
    - avg chamber pressure
    - peak chamber pressure
    - avg ISP
    - total impulse
    """


    if len(simRecord.points) < 2:
        print("Not enough simulation data to calculate summary.")
        return

    liquidMask = get_liquid_phase_mask(
        simRecord,
        liquidTankName
    )


    # overall averages
    avg_OF = time_average(simRecord,"mixture_ratio")
    avg_thrust = time_average(simRecord,"thrust_n")
    avg_chamber_pressure_pa = time_average(simRecord, "chamber_pressure_pa")
    avg_isp = time_average(simRecord,"isp_s")

    # liquid phase averages
    liquid_avg_OF = time_average(simRecord,"mixture_ratio", liquidMask)
    liquid_avg_thrust = time_average(simRecord,"thrust_n", liquidMask)
    liquid_avg_chamber_pressure_pa = time_average(simRecord, "chamber_pressure_pa", liquidMask)
    liquid_avg_isp = time_average(simRecord,"isp_s", liquidMask)

    # convert pressures to bar
    avg_chamber_pressure_bar = (
        avg_chamber_pressure_pa / 1e5
        if avg_chamber_pressure_pa is not None
        else None
    )

    liquid_avg_chamber_pressure_bar = (
        liquid_avg_chamber_pressure_pa / 1e5
        if liquid_avg_chamber_pressure_pa is not None
        else None
    )

    # overall stats
    impulse, burn_time = trap_integrate(simRecord,"thrust_n")
    liquid_impulse, liquid_burn_time = trap_integrate(simRecord,"thrust_n", liquidMask)
    peak_thrust = max(
        (
            point.engine.thrust_n
            for point in simRecord.points
            if point.engine is not None
            and point.engine.thrust_n is not None
            and math.isfinite(point.engine.thrust_n)
        ),
        default=None,
    )
    peak_chamber_pressure_pa = max(
        (
            point.engine.chamber_pressure_pa
            for point in simRecord.points
            if point.engine is not None
            and point.engine.chamber_pressure_pa is not None
            and math.isfinite(point.engine.chamber_pressure_pa)
        ),
        default=None,
    )
    peak_chamber_pressure_bar = (
        peak_chamber_pressure_pa / 1e5
        if peak_chamber_pressure_pa is not None
        else None
    )

    return SimSummary(
        burn_time_s=burn_time,
        liquid_burn_time_s=liquid_burn_time,

        avg_of=avg_OF,
        liquid_avg_of=liquid_avg_OF,

        avg_thrust_n=avg_thrust,
        liquid_avg_thrust_n=liquid_avg_thrust,
        peak_thrust_n=peak_thrust,

        avg_chamber_pressure_bar=avg_chamber_pressure_bar,
        liquid_avg_chamber_pressure_bar=liquid_avg_chamber_pressure_bar,
        peak_chamber_pressure_bar=peak_chamber_pressure_bar,

        avg_isp_s=avg_isp,
        liquid_avg_isp_s=liquid_avg_isp,

        total_impulse_ns=impulse,
        liquid_impulse_ns=liquid_impulse,
    )


def print_sim_summary(
    simRecord: SimRecord,
    liquidTankName: str,
) -> None:
    """
    Prints overall and liquid-phase engine performance stats.
    """

    summary = calculate_sim_summary(
        simRecord,
        liquidTankName,
    )

    if summary is None:
        print("Not enough simulation data to calculate summary.")
        return


    print()
    print("=== Engine Performance Summary ===")
    print(f"Burn time:             {summary.burn_time_s:10.3f} s")
    print(f"Liquid phase:          {summary.liquid_burn_time_s:10.3f} s")
    print()
    print("                         Overall       Liquid phase")
    print("----------------------------------------------------")
    print(
        f"Average O/F:          "
        f"{format_value(summary.avg_of):>10}       "
        f"{format_value(summary.liquid_avg_of):>10}"
    )
    print(
        f"Average thrust:       "
        f"{format_value(summary.avg_thrust_n):>10} N     "
        f"{format_value(summary.liquid_avg_thrust_n):>10} N"
    )
    print(
        f"Average chamber P:    "
        f"{format_value(summary.avg_chamber_pressure_bar):>10} bar   "
        f"{format_value(summary.liquid_avg_chamber_pressure_bar):>10} bar"
    )
    print(
        f"Average Isp:          "
        f"{format_value(summary.avg_isp_s):>10} s     "
        f"{format_value(summary.liquid_avg_isp_s):>10} s"
    )
    print()
    print(f"Total impulse:        {summary.total_impulse_ns:10.1f} Ns")
    print(f"Liquid impulse:       {summary.liquid_impulse_ns:10.1f} Ns")
    print(
        f"Peak thrust:          "
        f"{format_value(summary.peak_thrust_n):>10} N"
    )
    print(
        f"Peak chamber pressure:"
        f"{format_value(summary.peak_chamber_pressure_bar):>10} bar"
    )
    print("====================================================")



def plot_linkedin_summary(
    simRecord: SimRecord,
    oxidiserInjectorName: str,
    fuelInjectorName: str,
    file_path: str | Path | None = None,
    show: bool = True,
) -> None:
    """
    Creates a compact 2x2 simulation summary figure.

    Plots:
        - thrust
        - oxidiser tank mass breakdown
        - fuel tank mass breakdown
        - oxidiser and fuel injector mass flow
    """

    (
        timeSList,
        series,
        tankNames,
        injectorNames,
        _,
    ) = build_sim_series(simRecord)

    if len(series) == 0:
        print("No series to plot.")
        return

    firstPoint = simRecord.points[0]

    # --------------------------------------------------------
    # Identify fuel and oxidiser tanks
    # --------------------------------------------------------

    tankConfigs = {
        tankName: tankState.config
        for tankName, tankState in firstPoint.tanks.items()
    }

    oxidiserTankName, _ = get_single_config_by_role(
        tankConfigs,
        "oxidiser",
        "tank",
    )

    fuelTankName, _ = get_single_config_by_role(
        tankConfigs,
        "fuel",
        "tank",
    )

    # --------------------------------------------------------
    # Validate injector names
    # --------------------------------------------------------

    if oxidiserInjectorName not in injectorNames:
        raise ValueError(
            f"Oxidiser injector '{oxidiserInjectorName}' not found. "
            f"Available injectors: {injectorNames}"
        )

    if fuelInjectorName not in injectorNames:
        raise ValueError(
            f"Fuel injector '{fuelInjectorName}' not found. "
            f"Available injectors: {injectorNames}"
        )

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12, 8),
        sharex=True,
        constrained_layout=True,
    )

    thrustAx = axes[0, 0]
    oxidiserMassAx = axes[0, 1]
    fuelMassAx = axes[1, 0]
    mdotAx = axes[1, 1]

    # --------------------------------------------------------
    # Thrust
    # --------------------------------------------------------

    thrustAx.plot(
        timeSList,
        series["engine_thrust_n"],
    )
    thrustAx.set_title("Thrust")
    thrustAx.set_ylabel("Thrust (N)")
    thrustAx.set_ylim(bottom=0)

    # --------------------------------------------------------
    # Oxidiser tank mass
    # --------------------------------------------------------

    oxidiserTankState = firstPoint.tanks[oxidiserTankName]
    oxidiserPhaseModel = get_tank_phase_model(oxidiserTankState)

    plot_tank_mass(
        oxidiserMassAx,
        timeSList,
        series,
        oxidiserTankName,
        oxidiserPhaseModel,
    )

    oxidiserMassAx.set_title("N₂O Tank Mass")
    oxidiserMassAx.set_ylabel("Mass (kg)")

    # --------------------------------------------------------
    # Fuel tank mass
    # --------------------------------------------------------

    fuelTankState = firstPoint.tanks[fuelTankName]
    fuelPhaseModel = get_tank_phase_model(fuelTankState)

    plot_tank_mass(
        fuelMassAx,
        timeSList,
        series,
        fuelTankName,
        fuelPhaseModel,
    )

    fuelMassAx.set_title("Ethanol Tank Mass")
    fuelMassAx.set_ylabel("Mass (kg)")

    # --------------------------------------------------------
    # Propellant mass flow
    # --------------------------------------------------------

    mdotAx.plot(
        timeSList,
        series[f"{oxidiserInjectorName}_mdot_kg_s"],
        label="N₂O",
    )
    mdotAx.plot(
        timeSList,
        series[f"{fuelInjectorName}_mdot_kg_s"],
        label="Ethanol",
    )

    mdotAx.set_title("Propellant Mass Flow")
    mdotAx.set_ylabel("Mass flow (kg/s)")
    mdotAx.set_ylim(bottom=0)
    mdotAx.legend()

    # --------------------------------------------------------
    # Shared formatting
    # --------------------------------------------------------

    for ax in axes.flat:
        ax.grid(alpha=0.25)

    fuelMassAx.set_xlabel("Time (s)")
    mdotAx.set_xlabel("Time (s)")

    # --------------------------------------------------------
    # Save / display
    # --------------------------------------------------------

    if file_path is not None:
        fig.savefig(
            file_path,
            dpi=300,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    plt.close(fig)