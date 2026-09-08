from __future__ import annotations


from dataclasses import dataclass

from injector import InjectorConfig
from tanks import TankState, TankConfig
from regulator import RegulatorConfig
from constants import *


# ------------------------------------------------------------------
# Phase determination
#
# A saturated tank is two-phase whenever its bulk density sits above the
# saturated-vapour density at the tank temperature. This is re-evaluated from
# (mass, volume, temperature) every timestep and is never latched, so a tank
# that re-condenses after an apparent dryout is picked back up correctly.


@dataclass
class SaturatedSplit:
    is_two_phase: bool
    liquid_mass_kg: float
    vapour_mass_kg: float
    vapour_mass_fraction: float     # quality by mass, 0 = all liquid, 1 = all vapour
    saturation_properties: dict | None


def get_saturated_split(
    tank_config: TankConfig,
    tank_state: TankState
) -> SaturatedSplit:
    """
    Resolve the equilibrium liquid/vapour split of a saturated tank from its
    current mass, volume and temperature.
    """

    if tank_state.total_mass_kg is None or tank_state.temperature_k is None:
        raise ValueError("Tank state is missing total mass or temperature.")

    total_mass_kg = tank_state.total_mass_kg

    if total_mass_kg <= 0.0:
        return SaturatedSplit(False, 0.0, 0.0, 1.0, None)

    try:
        saturation_properties = tank_config.fluid.get_saturation_properties_from_temp(
            tank_state.temperature_k
        )
    except ValueError:
        # above the critical point or below the triple point - not saturated
        return SaturatedSplit(False, 0.0, total_mass_kg, 1.0, None)

    vf = saturation_properties["vf"]
    vg = saturation_properties["vg"]

    bulk_specific_volume_m3_kg = tank_config.tank_volume_m3 / total_mass_kg

    # bulk specific volume at or beyond the saturated-vapour line means there is
    # no liquid left and the fluid is superheated vapour
    if bulk_specific_volume_m3_kg >= vg:
        return SaturatedSplit(False, 0.0, total_mass_kg, 1.0, saturation_properties)

    vapour_mass_kg = (tank_config.tank_volume_m3 - total_mass_kg * vf) / (vg - vf)

    # clamp against round-off near the phase boundaries
    vapour_mass_kg = min(max(vapour_mass_kg, 0.0), total_mass_kg)
    liquid_mass_kg = total_mass_kg - vapour_mass_kg

    return SaturatedSplit(
        is_two_phase=True,
        liquid_mass_kg=liquid_mass_kg,
        vapour_mass_kg=vapour_mass_kg,
        vapour_mass_fraction=vapour_mass_kg / total_mass_kg,
        saturation_properties=saturation_properties,
    )


# Width of the blend band, in liquid volume fraction, over which the injector
# inlet transitions from pure liquid to pure vapour. A real tank does not switch
# instantaneously - as the liquid level reaches the outlet the injector ingests a
# rising vapour fraction. Blending over a narrow band also keeps the chamber
# pressure solve free of a step discontinuity at dryout.
LIQUID_HANDOVER_VOLUME_FRACTION = 0.02


def get_liquid_feed_weight(
    tank_config: TankConfig,
    split: SaturatedSplit
) -> float:
    """
    Fraction of the injector inlet flow that is liquid, 1.0 while the outlet is
    comfortably submerged and falling to 0.0 as the last liquid is drawn down.
    """

    if not split.is_two_phase:
        return 0.0

    if split.saturation_properties is None:
        return 0.0

    liquid_volume_m3 = split.liquid_mass_kg * split.saturation_properties["vf"]
    liquid_volume_fraction = liquid_volume_m3 / tank_config.tank_volume_m3

    if liquid_volume_fraction >= LIQUID_HANDOVER_VOLUME_FRACTION:
        return 1.0

    return liquid_volume_fraction / LIQUID_HANDOVER_VOLUME_FRACTION


def get_blowdown_injector_mdot_kg_s(
    tank_config: TankConfig,
    injector_config: InjectorConfig,
    downstream_pressure_pa: float,
    tank_state: TankState | None = None
) -> float:
    """
    Injector mass flow for a self-pressurised tank at its current state.

    Pure - advances nothing and mutates nothing, so the chamber pressure solver
    can evaluate it repeatedly at trial downstream pressures.
    """

    state = tank_state if tank_state is not None else tank_config.state

    if state is None:
        raise ValueError("Tank state is not initialised.")
    if state.total_mass_kg is None or state.total_mass_kg <= 0.0:
        return 0.0

    split = get_saturated_split(tank_config, state)
    liquid_weight = get_liquid_feed_weight(tank_config, split)

    if liquid_weight >= 1.0:
        return injector_config.get_dyer_mdot_kg_s(
            tank_state=state,
            downstream_pressure_pa=downstream_pressure_pa
        )

    gas_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
        tank_state=state,
        downstream_pressure_pa=downstream_pressure_pa
    )

    if liquid_weight <= 0.0:
        return gas_mdot_kg_s

    liquid_mdot_kg_s = injector_config.get_dyer_mdot_kg_s(
        tank_state=state,
        downstream_pressure_pa=downstream_pressure_pa
    )

    return liquid_weight * liquid_mdot_kg_s + (1.0 - liquid_weight) * gas_mdot_kg_s


def get_propellant_injector_mdot_kg_s(
    tank_config: TankConfig,
    injector_config: InjectorConfig,
    downstream_pressure_pa: float
) -> float:
    """
    Injector mass flow for any propellant tank at its current state, without
    advancing it. Used by the chamber pressure solver.
    """

    state = tank_config.state

    if state is None:
        raise ValueError("Tank state is not initialised.")

    if tank_config.phase_model == "self_pressurised":
        return get_blowdown_injector_mdot_kg_s(
            tank_config=tank_config,
            injector_config=injector_config,
            downstream_pressure_pa=downstream_pressure_pa
        )

    if tank_config.phase_model == "pressurised_liquid":

        if state.liquid_mass_kg is None or state.liquid_mass_kg <= 0.0:
            return 0.0

        return injector_config.get_liquid_mdot_kg_s(
            tank_state=state,
            downstream_pressure_pa=downstream_pressure_pa
        )

    if tank_config.phase_model == "single_phase":
        return injector_config.get_gas_mdot_kg_s(
            tank_state=state,
            downstream_pressure_pa=downstream_pressure_pa
        )

    raise NotImplementedError(
        f"Unsupported phase model '{tank_config.phase_model}' for tank {tank_config.name}"
    )


def gas_advance_timestep(
    tank_config: TankConfig,
    injector_config: InjectorConfig,
    dt_s: float,
    downstream_pressure_pa: float = ATMOSPHERE_PRESSURE_PA
) -> tuple[TankState, float]:
    """
    Advances the state of a gas-only tank by one timestep of length dt_s. 
    """


    if tank_config.state is None:
        raise ValueError("Tank state is not initialised.")
    if tank_config.state.pressure_pa is None or tank_config.state.temperature_k is None:
        raise ValueError("Tank state is missing pressure or temperature.")
    if tank_config.state.total_mass_kg is None or tank_config.state.total_internal_energy_j is None:
        raise ValueError("Tank state is missing total mass or total internal energy.")


    # advancing the gas timestep is much simplier than the self-pressurised tank, as there's only 1 phase the fluid can be in, so we can always safely default to gas injection

    # find current tank pressure & temp
    tank_pressure_pa = tank_config.state.pressure_pa
    tank_temperature_k = tank_config.state.temperature_k


    # find the mdot through the injector
    injector_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
        tank_state=tank_config.state,
        downstream_pressure_pa=downstream_pressure_pa
    )

    # total mass leaving the tank this timestep
    mass_out_kg = injector_mdot_kg_s*dt_s

    # we can find the total energy left in the tank after this timestep by finding how much energy leaves the system via the leaving mass
    # for a real gas-only state we aren't saturated, so need to use actual gas enthalpy
    gas_density_kg_m3 = tank_config.state.total_mass_kg / tank_config.tank_volume_m3

    gas_enthalpy_j_kg = tank_config.fluid.props_si(
        "H", "D", gas_density_kg_m3, "T", tank_temperature_k
    )

    energy_out_j = mass_out_kg * gas_enthalpy_j_kg

    # find new total mass and energy
    new_total_mass_kg = tank_config.state.total_mass_kg - mass_out_kg
    new_total_energy_j = tank_config.state.total_internal_energy_j - energy_out_j

    # use this to determine the new state of the tank
    new_tank_state = tank_config.state_from_mass_and_energy(
        total_mass_kg=new_total_mass_kg,
        total_internal_energy_j=new_total_energy_j,
        previous_state=tank_config.state,
        phase_override="single_phase"
    )


    return new_tank_state, injector_mdot_kg_s



def blowdown_advance_timestep(
    tank_config: TankConfig,
    injector_config: InjectorConfig,
    dt_s: float,
    downstream_pressure_pa: float = ATMOSPHERE_PRESSURE_PA
) -> tuple[TankState, float]:
    """
    Advance a self-pressurised tank by one timestep.

    The phase of the tank is resolved every timestep from its bulk density
    against the saturated-vapour density at the tank temperature. Nothing is
    latched: a tank that re-condenses after an apparent dryout is returned to the
    two-phase model automatically.

    Mass and energy always leave on the same assumption - the liquid/vapour split
    of the outflow follows the same weighting used to blend the injector models,
    so the enthalpy debited matches the fluid actually expelled.
    """

    current_state = tank_config.state

    if current_state is None:
        raise ValueError("Tank state is not initialised.")
    if current_state.pressure_pa is None or current_state.temperature_k is None:
        raise ValueError("Tank state is missing pressure or temperature.")
    if current_state.total_mass_kg is None or current_state.total_internal_energy_j is None:
        raise ValueError("Tank state is missing total mass or total internal energy.")

    if current_state.total_mass_kg <= 0.0:
        return current_state, 0.0

    split = get_saturated_split(tank_config, current_state)
    liquid_weight = get_liquid_feed_weight(tank_config, split)

    mdot_kg_s = get_blowdown_injector_mdot_kg_s(
        tank_config=tank_config,
        injector_config=injector_config,
        downstream_pressure_pa=downstream_pressure_pa,
        tank_state=current_state
    )

    if mdot_kg_s <= 0.0:
        return current_state, 0.0

    mass_out_kg = mdot_kg_s * dt_s

    # never remove more mass than the tank holds
    if mass_out_kg > current_state.total_mass_kg:
        mass_out_kg = current_state.total_mass_kg
        mdot_kg_s = mass_out_kg / dt_s

    # ------------------------------------------------------------------
    # energy carried out by the departing mass
    #
    # The outflow is liquid_weight liquid and (1 - liquid_weight) vapour, matching
    # how the injector models were blended, so mass and energy stay consistent.

    if split.is_two_phase and split.saturation_properties is not None:

        saturation_properties = split.saturation_properties

        outlet_enthalpy_j_kg = (
            liquid_weight * saturation_properties["hf"]
            + (1.0 - liquid_weight) * saturation_properties["hg"]
        )

    else:
        # superheated vapour - use the real gas enthalpy at the bulk state
        bulk_density_kg_m3 = current_state.total_mass_kg / tank_config.tank_volume_m3

        outlet_enthalpy_j_kg = tank_config.fluid.get_fluid_enthalpy_from_density_temperature(
            D=bulk_density_kg_m3,
            T=current_state.temperature_k
        )

    energy_out_j = mass_out_kg * outlet_enthalpy_j_kg

    new_total_mass_kg = current_state.total_mass_kg - mass_out_kg
    new_total_internal_energy_j = current_state.total_internal_energy_j - energy_out_j

    # ------------------------------------------------------------------
    # reconstruct the new state, choosing the phase model from the result rather
    # than from a persisted flag

    new_tank_state = reconstruct_blowdown_state(
        tank_config=tank_config,
        total_mass_kg=new_total_mass_kg,
        total_internal_energy_j=new_total_internal_energy_j,
        previous_state=current_state
    )

    return new_tank_state, mdot_kg_s


def reconstruct_blowdown_state(
    tank_config: TankConfig,
    total_mass_kg: float,
    total_internal_energy_j: float,
    previous_state: TankState
) -> TankState:
    """
    Rebuild a self-pressurised tank state from mass and energy, selecting the
    two-phase or single-phase solver based on which one the result is consistent
    with. Falls back to the single-phase solver if the saturated solve cannot
    represent the state.
    """

    if total_mass_kg <= 0.0:
        raise ValueError(f"Tank {tank_config.name} ran out of mass.")

    bulk_specific_volume_m3_kg = tank_config.tank_volume_m3 / total_mass_kg

    # a state can only be two-phase if its bulk specific volume falls between the
    # saturated liquid and saturated vapour lines somewhere on the dome
    try:
        saturated_state = tank_config.state_from_mass_and_energy(
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            previous_state=previous_state,
            phase_override="self_pressurised"
        )

        if (
            saturated_state.temperature_k is not None
            and saturated_state.liquid_mass_kg is not None
            and saturated_state.liquid_mass_kg > 0.0
        ):
            saturation_properties = tank_config.fluid.get_saturation_properties_from_temp(
                saturated_state.temperature_k
            )

            if bulk_specific_volume_m3_kg < saturation_properties["vg"]:
                return saturated_state

    except (ValueError, RuntimeError):
        pass

    # no valid saturated state - the tank is superheated vapour
    return tank_config.state_from_mass_and_energy(
        total_mass_kg=total_mass_kg,
        total_internal_energy_j=total_internal_energy_j,
        previous_state=previous_state,
        phase_override="single_phase"
    )


def advance_pressurant_tank(
        tank_config: TankConfig,
        dt_s: float,
        pressurant_mdot_kg_s: float
) -> TankState:
    

    if tank_config.state is None:
        raise ValueError("tank_config is missing state")
    if tank_config.state.total_mass_kg is None:
        raise ValueError("tank_state is missing total_mass_kg")
    if tank_config.state.total_internal_energy_j is None:
        raise ValueError("tank_state is missing total_internal_energy_j")
    if tank_config.state.temperature_k is None:
        raise ValueError("tank_state is missing temperature_k")
    
    


    fluid_density_kg_m3 = tank_config.state.total_mass_kg / tank_config.tank_volume_m3

    pressurant_enthalpy_j_kg = tank_config.fluid.get_fluid_enthalpy_from_density_temperature(
        D=fluid_density_kg_m3,
        T=tank_config.state.temperature_k
    )

    mass_out_kg = pressurant_mdot_kg_s * dt_s
    energy_out_j = pressurant_mdot_kg_s * dt_s * pressurant_enthalpy_j_kg

    new_total_mass_kg = tank_config.state.total_mass_kg - mass_out_kg
    new_total_internal_energy_j = tank_config.state.total_internal_energy_j - energy_out_j



    new_pressurant_tank_state = tank_config.state_from_mass_and_energy(
        total_mass_kg=new_total_mass_kg,
        total_internal_energy_j=new_total_internal_energy_j,
        previous_state=tank_config.state,
        phase_override="single_phase"
    )

    return new_pressurant_tank_state





def pressurised_liquid_advance_timestep(
        liquid_tank_config: TankConfig,
        injector_config: InjectorConfig,
        regulator_config: RegulatorConfig,
        pressurant_tank_config: TankConfig,
        dt_s: float,
        downstream_pressure_pa: float = ATMOSPHERE_PRESSURE_PA
) -> tuple[ TankState, float, float]:
    

    liquid_state = liquid_tank_config.state
    if liquid_state is None:
        raise ValueError("Liquid tank state is not initialised.")
    if liquid_state.pressure_pa is None or liquid_state.temperature_k is None:
        raise ValueError("Liquid tank state is missing pressure or temperature.")
    if liquid_state.total_mass_kg is None:
        raise ValueError("Liquid tank state is missing total mass.")
    if liquid_state.total_internal_energy_j is None:
        raise ValueError("Liquid tank state is missing total internal energy.")
    if liquid_state.pressurant_gas_mass_kg is None:
        raise ValueError("Liquid tank state is missing pressurant gas mass.")
    if liquid_state.liquid_mass_kg is None:
        raise ValueError("Liquid tank state is missing liquid mass.")
    
    pressurant_state = pressurant_tank_config.state
    if pressurant_state is None:
        raise ValueError("Pressurant tank state is not initialised.")
    if pressurant_state.pressure_pa is None or pressurant_state.temperature_k is None:
        raise ValueError("Pressurant tank state is missing pressure or temperature.")
    if pressurant_state.total_mass_kg is None:
        raise ValueError("Pressurant tank state is missing total mass.")
    if pressurant_state.total_internal_energy_j is None:
        raise ValueError("Pressurant tank state is missing total internal energy.")
    
    if liquid_tank_config.pressurant_fluid is None:
        raise ValueError("Liquid tank config is missing a defined pressurant fluid.")
    
    #-------------------------
    # 1. remove liquid from the liquid tank via the injector

    liquid_pressure_pa = liquid_state.pressure_pa
    liquid_temperature_k = liquid_state.temperature_k
    liquid_mass_kg = liquid_state.liquid_mass_kg
    liquid_pressurant_mass_kg = liquid_state.pressurant_gas_mass_kg
    liquid_tank_total_internal_energy_j = liquid_state.total_internal_energy_j

    liquid_density_kg_m3 = liquid_tank_config.fluid.get_fluid_density_from_pressure_temperature(
        P=liquid_pressure_pa,
        T=liquid_temperature_k
    )

    delta_p_pa = liquid_pressure_pa - downstream_pressure_pa
    if delta_p_pa <= 0.0 or liquid_mass_kg <= 0.0:
        # no flow through the injector
        liquid_mdot_kg_s = 0.0

    else:
        liquid_mdot_kg_s = injector_config.get_liquid_mdot_kg_s(
            tank_state=liquid_state,
            downstream_pressure_pa=downstream_pressure_pa
        )
    
    liquid_mass_out_kg = liquid_mdot_kg_s * dt_s

    # clamp so we never remove more liquid mass than is actually in the tank
    if liquid_mass_out_kg > liquid_mass_kg:
        liquid_mass_out_kg = liquid_mass_kg
        liquid_mdot_kg_s = liquid_mass_out_kg / dt_s
    
    # for energy out, we can use the enthalpy of the liquid at current tank pressure and temperature, as the flow is liquid at the injector inlet
    liquid_enthalpy_j_kg = liquid_tank_config.fluid.get_fluid_enthalpy_from_pressure_temperature(
        P=liquid_pressure_pa,
        T=liquid_temperature_k
    )
    liquid_energy_out_j = liquid_mass_out_kg * liquid_enthalpy_j_kg

    new_liquid_mass_kg = liquid_mass_kg - liquid_mass_out_kg

    liquid_dryout = new_liquid_mass_kg <= DRYOUT_TOLERANCE_KG
    if liquid_dryout:
        new_liquid_mass_kg = 0.0

    
    #-------------------------
    # 2. add persurant through the regulator

    pressurant_mdot_kg_s = regulator_config.get_regulator_mdot_kg_s(
        liquid_tank_state=liquid_state,
        pressurant_tank_state=pressurant_state,
        new_liquid_mass_kg=new_liquid_mass_kg,
        dt_s=dt_s
    )

    pressurant_mass_in_kg = pressurant_mdot_kg_s * dt_s

    # clamp againt the available mass in the pressurant tank
    if pressurant_mass_in_kg > pressurant_state.total_mass_kg:
        pressurant_mass_in_kg = pressurant_state.total_mass_kg
        pressurant_mdot_kg_s = pressurant_mass_in_kg / dt_s
    
    pressurant_tank_density_kg_m3 = pressurant_state.total_mass_kg / pressurant_tank_config.tank_volume_m3

    pressurant_inlet_enthalpy_j_kg = liquid_tank_config.pressurant_fluid.get_fluid_enthalpy_from_density_temperature(
        D=pressurant_tank_density_kg_m3,
        T=pressurant_state.temperature_k
    )

    pressurant_energy_in_j = pressurant_mass_in_kg * pressurant_inlet_enthalpy_j_kg

    new_liquid_pressurant_mass_kg = liquid_pressurant_mass_kg + pressurant_mass_in_kg

    #-------------------------
    # 3. update total internal energies

    new_liquid_tank_total_internal_energy_j = liquid_tank_total_internal_energy_j - liquid_energy_out_j + pressurant_energy_in_j

    # new_pressurant_tank_total_mass = pressurant_state.total_mass_kg - pressurant_mass_in_kg
    # new_pressurant_tank_total_internal_energy_j = pressurant_state.total_internal_energy_j - pressurant_energy_in_j


    # if new_pressurant_tank_total_mass <= 0.0:
    #     raise ValueError("Pressurant tank mass went negative, this should have been clamped against.")
    

    # -------------------------
    # 4. reconstruct tank states

    new_liquid_tank_state = liquid_tank_config.state_from_mass_and_energy(
        liquid_mass_kg=new_liquid_mass_kg,
        pressurant_gas_mass_kg=new_liquid_pressurant_mass_kg,
        total_internal_energy_j=new_liquid_tank_total_internal_energy_j,
        previous_state=liquid_state,
        phase_override="pressurised_liquid",
    )

    # new_pressurant_tank_state = pressurant_tank_config.state_from_mass_and_energy(
    #     total_mass_kg=new_pressurant_tank_total_mass,
    #     total_internal_energy_j=new_pressurant_tank_total_internal_energy_j,
    #     previous_state=pressurant_state,
    #     phase_override="single_phase"
    # )

    return (
        new_liquid_tank_state,
        liquid_mdot_kg_s,
        pressurant_mdot_kg_s
    )




# ------------
# General wrapper for advancing tank timesteps

@dataclass
class PropellantAdvanceResult:
    tank_update: TankState
    injector_mdot_kg_s: float | None
    regulator_mdot_kg_s: float | None


def advance_propellant_tank(
    tank_config: TankConfig,
    injector_config: InjectorConfig,
    dt_s: float,
    downstream_pressure_pa: float,
    pressurant_tank_config: TankConfig | None = None,
    regulator_config: RegulatorConfig | None = None
) -> PropellantAdvanceResult:
    

    if tank_config.phase_model == "self_pressurised":
        
        new_tank_state, injector_mdot_kg_s = blowdown_advance_timestep(
            tank_config=tank_config,
            injector_config=injector_config,
            dt_s=dt_s,
            downstream_pressure_pa=downstream_pressure_pa
        )


        return PropellantAdvanceResult(
            tank_update=new_tank_state,
            injector_mdot_kg_s=injector_mdot_kg_s,
            regulator_mdot_kg_s=None
            )
    

    if tank_config.phase_model == "pressurised_liquid":

        if (
            pressurant_tank_config is None
            or regulator_config is None
        ):
            raise ValueError(
                f"Tank '{tank_config.name}' uses phase_model = 'pressurised_liquid', so it requires a pressurant tank and a regulator."
            )
        

        (
            new_tank_state,
            injector_mdot_kg_s,
            regulator_mdot_kg_s
        ) = pressurised_liquid_advance_timestep(
            liquid_tank_config=tank_config,
            injector_config=injector_config,
            regulator_config=regulator_config,
            pressurant_tank_config=pressurant_tank_config,
            dt_s=dt_s,
            downstream_pressure_pa=downstream_pressure_pa
        )

        return PropellantAdvanceResult(
            tank_update=new_tank_state,
            injector_mdot_kg_s=injector_mdot_kg_s,
            regulator_mdot_kg_s=regulator_mdot_kg_s
        )
    


    raise NotImplementedError(
        f"Unsupported phase model '{tank_config.phase_model}' for tank {tank_config.name}"
    )
