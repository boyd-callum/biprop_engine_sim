
from dataclasses import dataclass
from fluid import Fluid
from tanks import TankState


@dataclass
class RegulatorConfig:
    name: str
    set_pressure_pa: float

    def get_regulator_mdot_kg_s(
            self,
            liquid_tank_state: TankState,
            pressurant_tank_state: TankState,
            new_liquid_mass_kg: float,
            dt_s: float
            ) -> float:

        """
        Find the mass flow rate of pressurant gas required to bring a liquid tank up to the regulator set pressure over the current timestep
        """

        if liquid_tank_state.config.pressurant_fluid is None:
            raise ValueError(f"Liquid tank ({liquid_tank_state.config.name}) requires a pressurant fluid.")
        if liquid_tank_state.config.pressurant_fluid != pressurant_tank_state.config.fluid:
            raise ValueError(f"Liquid tank ({liquid_tank_state.config.name}) pressurant fluid ({liquid_tank_state.config.pressurant_fluid}) is not the same as Pressurant tank ({pressurant_tank_state.config.name}) fluid ({pressurant_tank_state.config.fluid}).")
        

        # check inputs 
        if liquid_tank_state.ullage_volume_m3 is None:
            raise ValueError("Liquid tank state is missing ullage volume.")
        if liquid_tank_state.pressurant_gas_mass_kg is None:
            raise ValueError("Liquid tank state is missing pressurant gas mass.")
        if liquid_tank_state.temperature_k is None:
            raise ValueError("Liquid tank state is missing temperature.")
        

        if liquid_tank_state.config.pressurant_fluid is None:
            raise ValueError("Liquid tank config is missing a defined pressurant gas")

        pressurant_fluid = liquid_tank_state.config.pressurant_fluid


        downstream_temperature_k = liquid_tank_state.temperature_k
        current_pressurant_mass_kg = liquid_tank_state.pressurant_gas_mass_kg


        # estimate the liquid volume in the liquid tank at the regulator set pressure
        liquid_density_kg_m3 = liquid_tank_state.config.fluid.get_fluid_density_from_pressure_temperature(
            P=self.set_pressure_pa,
            T=downstream_temperature_k
        )
        liquid_volume_m3 = new_liquid_mass_kg / liquid_density_kg_m3
        ullage_volume_m3 = liquid_tank_state.config.tank_volume_m3 - liquid_volume_m3

        if ullage_volume_m3 <= 0.0:
            raise ValueError(
                f"Computed ullage volume in {liquid_tank_state.config.name} must remain positive."
            )


        # find required density and mass to meet target pressure
        required_pressurant_density_kg_m3 = pressurant_fluid.get_fluid_density_from_pressure_temperature(
            P=self.set_pressure_pa,
            T=downstream_temperature_k
        )
        required_pressurant_mass_kg = required_pressurant_density_kg_m3 * ullage_volume_m3


        # find difference between needed and current mass
        delta_mass_kg = required_pressurant_mass_kg - current_pressurant_mass_kg

        if delta_mass_kg <= 0.0:
            return 0.0
        
        return delta_mass_kg / dt_s
