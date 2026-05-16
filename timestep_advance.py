

from dataclasses import dataclass

from injector import InjectorConfig
from tanks import TankState, TankConfig
from regulator import RegulatorConfig
from constants import *



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

    current_state = tank_config.state

    if current_state is None:
        raise ValueError("Tank state is not initialised.")
    if current_state.pressure_pa is None or current_state.temperature_k is None:
        raise ValueError("Tank state is missing pressure or temperature.")
    if current_state.total_mass_kg is None or current_state.total_internal_energy_j is None:
        raise ValueError("Tank state is missing total mass or total internal energy.")
    if current_state.liquid_mass_kg is None:
        raise ValueError("Tank state is missing liquid mass.")

    

    tank_temperature_k = current_state.temperature_k
    liquid_mass_kg = current_state.liquid_mass_kg

    def advance_gas_from_state(
        start_state: TankState,
        gas_dt_s: float
    ) -> tuple[TankState, float, float, float]:
        """
        Advances a gas-only state by gas_dt_s.

        Returns:
            new_state, gas_mdot_kg_s, gas_mass_out_kg, gas_energy_out_j
        """

        if gas_dt_s <= 0.0:
            return start_state, 0.0, 0.0, 0.0

        if start_state.total_mass_kg is None:
            raise ValueError("Gas state is missing total mass.")
        if start_state.total_internal_energy_j is None:
            raise ValueError("Gas state is missing total internal energy.")
        if start_state.temperature_k is None:
            raise ValueError("Gas state is missing temperature.")
        


        if start_state.total_mass_kg <= 0.0:
            return start_state, 0.0, 0.0, 0.0

        gas_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
            tank_state=start_state,
            downstream_pressure_pa=downstream_pressure_pa
        )

        gas_mass_out_kg = gas_mdot_kg_s * gas_dt_s

        # Clamp so we never remove more mass than exists.
        if gas_mass_out_kg > start_state.total_mass_kg:
            gas_mass_out_kg = start_state.total_mass_kg
            gas_mdot_kg_s = gas_mass_out_kg / gas_dt_s

        gas_density_kg_m3 = start_state.total_mass_kg / tank_config.tank_volume_m3

        gas_enthalpy_j_kg = tank_config.fluid.props_si(
            "H",
            "D",
            gas_density_kg_m3,
            "T",
            start_state.temperature_k
        )

        gas_energy_out_j = gas_mass_out_kg * gas_enthalpy_j_kg

        new_total_mass_kg = start_state.total_mass_kg - gas_mass_out_kg
        new_total_internal_energy_j = (
            start_state.total_internal_energy_j
            - gas_energy_out_j
        )

        new_state = tank_config.state_from_mass_and_energy(
            total_mass_kg=new_total_mass_kg,
            total_internal_energy_j=new_total_internal_energy_j,
            previous_state=start_state,
            phase_override="single_phase"
        )

        return new_state, gas_mdot_kg_s, gas_mass_out_kg, gas_energy_out_j

    # ------------------------------------------------------------------
    # Case 1: tank is already gas-only at the start of the timestep

    if liquid_mass_kg <= DRYOUT_TOLERANCE_KG:

        new_tank_state, gas_mdot_kg_s, _, _ = advance_gas_from_state(
            start_state=current_state,
            gas_dt_s=dt_s
        )

        return new_tank_state, gas_mdot_kg_s

    # ------------------------------------------------------------------
    # Case 2: liquid exists at the start of the timestep

    saturation_properties = tank_config.fluid.get_saturation_properties_from_temp(
        tank_temperature_k
    )

    liquid_mdot_kg_s = injector_config.get_dyer_mdot_kg_s(
        tank_state=current_state,
        downstream_pressure_pa=downstream_pressure_pa
    )

    if liquid_mdot_kg_s <= 0.0:
        return current_state, 0.0

    full_step_liquid_mass_out_kg = liquid_mdot_kg_s * dt_s

    # ------------------------------------------------------------------
    # Case 2a: liquid remains for the whole timestep

    if liquid_mass_kg > full_step_liquid_mass_out_kg + DRYOUT_TOLERANCE_KG:

        mass_out_kg = full_step_liquid_mass_out_kg
        energy_out_j = mass_out_kg * saturation_properties["hf"]

        new_total_mass_kg = current_state.total_mass_kg - mass_out_kg
        new_total_internal_energy_j = (
            current_state.total_internal_energy_j
            - energy_out_j
        )

        new_tank_state = tank_config.state_from_mass_and_energy(
            total_mass_kg=new_total_mass_kg,
            total_internal_energy_j=new_total_internal_energy_j,
            previous_state=current_state,
            phase_override="self_pressurised"
        )

        return new_tank_state, liquid_mdot_kg_s

    # ------------------------------------------------------------------
    # Case 2b: dryout occurs inside this timestep

    liquid_time_s = liquid_mass_kg / liquid_mdot_kg_s
    liquid_time_s = max(0.0, min(liquid_time_s, dt_s))

    gas_time_s = dt_s - liquid_time_s

    # First remove exactly the remaining liquid.
    liquid_mass_out_kg = liquid_mass_kg
    liquid_energy_out_j = liquid_mass_out_kg * saturation_properties["hf"]

    dryout_total_mass_kg = current_state.total_mass_kg - liquid_mass_out_kg
    dryout_total_internal_energy_j = (
        current_state.total_internal_energy_j
        - liquid_energy_out_j
    )

    # Reconstruct the intermediate dryout state as gas-only.
    dryout_state = tank_config.state_from_mass_and_energy(
        total_mass_kg=dryout_total_mass_kg,
        total_internal_energy_j=dryout_total_internal_energy_j,
        previous_state=current_state,
        phase_override="single_phase"
    )

    # Since this is explicitly the dryout point, do not allow a tiny
    # reconstructed liquid mass to carry through this timestep.
    dryout_state.liquid_mass_kg = 0.0

    # Then use gas discharge for the remaining part of the timestep.
    (
        final_state,
        gas_mdot_kg_s,
        gas_mass_out_kg,
        gas_energy_out_j
    ) = advance_gas_from_state(
        start_state=dryout_state,
        gas_dt_s=gas_time_s
    )

    final_state.liquid_mass_kg = 0.0

    total_mass_out_kg = liquid_mass_out_kg + gas_mass_out_kg

    # Return timestep-averaged injector mdot so the rest of the sim does not
    # need to know that this timestep was internally split.
    average_injector_mdot_kg_s = total_mass_out_kg / dt_s

    return final_state, average_injector_mdot_kg_s
    




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
