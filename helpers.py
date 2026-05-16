from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Mapping, overload, Callable

if TYPE_CHECKING:
    from simulation import SimRecord
    from tanks import TankConfig, SourceRole
    from injector import InjectorConfig
    from regulator import RegulatorConfig, PropellantRole


def bisection_search(
        residual_func: Callable, 
        bounds: list, 
        max_iterations: int = 100,
        tolerance: float = 1e-3
        ) -> float:

    low = bounds[0]
    high = bounds[1]

    low_residual = residual_func(low)
    high_residual = residual_func(high)

    if low_residual * high_residual > 0.0:
        raise ValueError(
            f"Bisection bounds do not bracket a root: "
            f"f({low})={low_residual}, f({high})={high_residual}"
        )

    iteration = 0

    while iteration < max_iterations:
        
        mid = 0.5 * (low + high)

        mid_residual = residual_func(mid)
        
        # stop once residual is close enough to zero
        if abs(mid_residual) < tolerance:
            break

        # keep the half-interval that still has the sign change
        if low_residual * mid_residual <= 0.0:
            high = mid
            high_residual = mid_residual
        else:
            low = mid
            low_residual = mid_residual
        
        iteration +=1

        if iteration > 100:
            raise RuntimeError(f"Root search failed to converge after 100 iterations. Final residual of {mid_residual}")

    return 0.5 * (low + high)





ConfigTypes = Literal["tank", "injector", "regulator"]

# we use these weird overload typing hints so that the linter knows what the returned type will be. this is because the returned type depends on the inputed value for the config type

# eg. if input is "tank", then output is TankConfig

@overload
def get_single_config_by_role(
    configs: Mapping[str, TankConfig],
    role: SourceRole,
    config_type_name: Literal["tank"],
) -> tuple[str, TankConfig]:
    ...

@overload
def get_single_config_by_role(
    configs: Mapping[str, InjectorConfig],
    role: PropellantRole,
    config_type_name: Literal["injector"],
) -> tuple[str, InjectorConfig]:
    ...


@overload
def get_single_config_by_role(
    configs: Mapping[str, RegulatorConfig],
    role: PropellantRole,
    config_type_name: Literal["regulator"],
) -> tuple[str, RegulatorConfig]:
    ...


def get_single_config_by_role(
    configs: (
        Mapping[str, TankConfig]
        | Mapping[str, InjectorConfig]
        | Mapping[str, RegulatorConfig]
    ),
    role: SourceRole,
    config_type_name: ConfigTypes,
) -> tuple[str, TankConfig | InjectorConfig | RegulatorConfig]:
    """
    Finds one config object with the specified role.

    Used when a sim requires a single component for a given role, such as one fuel tank, one oxidiser tank, or one fuel injector. Searches through a dictionary of config objects and returns matching dictionary key and congig object.
    """
    
    if (config_type_name == "injector" or config_type_name == "regulator") and role == "pressurant":
        raise ValueError(f"{config_type_name} config cannot have role 'pressurant'.")

    matching_items = [
        (config_id, config)
        for config_id, config in configs.items()
        if config.role == role
    ]

    if len(matching_items) == 0:
        raise ValueError(f"No {config_type_name} config found with role '{role}'.")

    if len(matching_items) > 1:
        matching_ids = [config_id for config_id, _ in matching_items]
        raise ValueError(
            f"Multiple {config_type_name}configs found with role '{role}': "
            f"{matching_ids}"
        )

    return matching_items[0]



def format_bar(value_pa: float | None) -> str:
    if value_pa is None:
        return "N/A"
    return f"{value_pa / 1e5:.2f}"


def format_force(value_n: float | None) -> str:
    if value_n is None:
        return "N/A"
    return f"{value_n:.1f}"




def get_tank_debug_string(tank_id: str, tank_config: TankConfig) -> str:
    """
    Builds a compact tank state string for exception context.
    """

    state = tank_config.state

    if state is None:
        return f"{tank_id}: state=None"

    return (
        f"{tank_id}: "
        f"name={tank_config.name}, "
        f"role={tank_config.role}, "
        f"phase_model={tank_config.phase_model}, "
        f"P={format_bar(state.pressure_pa)}, "
        f"T={state.temperature_k}, "
        f"liquid_mass={state.liquid_mass_kg}, "
        f"vapour_mass={state.vapour_mass_kg}, "
        f"pressurant_mass={state.pressurant_gas_mass_kg}, "
        f"total_mass={state.total_mass_kg}"
    )