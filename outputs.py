from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt

import csv

from simulation import SimRecord


def safe_divide(value: float | None, divisor: float) -> float | None:
    """
    Safely divides a value by a divisor, preserving None values.
    """

    if value is None:
        return None

    return value / divisor


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


def get_ylabel(seriesKey: str) -> str:
    """
    Returns a suitable y-axis label for a plotted series key.
    """

    if "pressure" in seriesKey:
        return "Pressure (bar)"
    if "temperature" in seriesKey:
        return "Temperature (K)"
    if "mass" in seriesKey:
        return "Mass (kg)"
    if "energy" in seriesKey:
        return "Energy (kJ)"
    if "mdot" in seriesKey:
        return "Mass flow rate (kg/s)"
    if "thrust" in seriesKey:
        return "Thrust (N)"
    if "isp" in seriesKey:
        return "Specific impulse (s)"
    if "mixture_ratio" in seriesKey:
        return "O/F ratio"
    if "cstar" in seriesKey:
        return "c* (m/s)"
    if "cf" in seriesKey:
        return "Thrust coefficient"

    return "Value"


def should_force_positive_y(seriesKey: str) -> bool:
    """
    Returns True if the plotted quantity should normally be non-negative.
    """

    positiveKeywords = [
        "pressure",
        "mass",
        "mdot",
        "thrust",
        "isp",
        "cstar",
        "cf",
    ]

    return any(keyword in seriesKey for keyword in positiveKeywords)


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


def get_combined_mass_plot_keys(
    simRecord: SimRecord,
    tankNames: list[str],
) -> list[str]:
    """
    Returns combined mass plot keys for tanks where multiple mass components
    are useful to view together.
    """

    if len(simRecord.points) == 0:
        return []

    firstPoint = simRecord.points[0]

    massModels = {"self_pressurised", "pressurised_liquid"}
    massPlotKeys: list[str] = []

    for tankName in tankNames:
        tankState = firstPoint.tanks.get(tankName) if firstPoint.tanks else None

        if tankState is None:
            continue

        phaseModel = getattr(tankState.config, "phase_model", None)

        if phaseModel in massModels:
            massPlotKeys.append(f"{tankName}_all_masses")

    return massPlotKeys


def plot_sim_record(
    simRecord: SimRecord,
    file_path: str | Path | None = None,
    cols: int = 3,
    show: bool = True,
) -> None:
    """
    Plots all available tank, injector, regulator, and engine data in a
    SimRecord.

    Args:
        simRecord:
            Completed simulation record.
        savePath:
            Optional file path for saving the generated figure.
        cols:
            Number of subplot columns.
        show:
            Whether to display the figure using ``plt.show()``.
    """

    timeSList, series, tankNames, _, _ = build_sim_series(simRecord)

    if len(series) == 0:
        print("No series to plot.")
        return

    plotKeys = list(series.keys())
    plotKeys.extend(get_combined_mass_plot_keys(simRecord, tankNames))

    numPlots = len(plotKeys)
    rows = (numPlots + cols - 1) // cols

    fig, axs = plt.subplots(
        rows,
        cols,
        figsize=(4 * cols, 3 * rows),
    )

    axsFlat = axs.flatten() if hasattr(axs, "flatten") else [axs]

    lastIndex = -1

    for lastIndex, plotKey in enumerate(plotKeys):
        ax = axsFlat[lastIndex]

        if plotKey.endswith("_all_masses"):
            tankName = plotKey.removesuffix("_all_masses")

            ax.plot(
                timeSList,
                series[f"{tankName}_liquid_mass_kg"],
                label="liquid mass",
            )
            ax.plot(
                timeSList,
                series[f"{tankName}_vapour_mass_kg"],
                label="vapour mass",
            )
            ax.plot(
                timeSList,
                series[f"{tankName}_pressurant_gas_mass_kg"],
                label="pressurant mass",
            )
            ax.plot(
                timeSList,
                series[f"{tankName}_total_mass_kg"],
                label="total mass",
            )

            ax.set_title(f"{tankName} all masses")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Mass (kg)")
            ax.set_ylim(bottom=0)
            ax.legend()
            continue

        ax.plot(timeSList, series[plotKey])
        ax.set_title(plotKey.replace("_", " "))
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(get_ylabel(plotKey))

        if should_force_positive_y(plotKey):
            ax.set_ylim(bottom=0)

    for axisIndex in range(lastIndex + 1, len(axsFlat)):
        axsFlat[axisIndex].axis("off")

    plt.tight_layout()

    if file_path is not None:
        fig.savefig(file_path, dpi=300)

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