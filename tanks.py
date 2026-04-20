

from dataclasses import dataclass
from typing import Literal

from injector import InjectorConfig
from fluid import Fluid

# -------------------------------
# Dataclasses for the tanks

SourceRole = Literal["fuel", "oxidiser", "pressurant"]
PhaseModel = Literal["liquid", "gas", "self_pressurised", "unknown"]
TankInitMode = Literal[
    "pressure_mass",         # good for self-pressurised N2O
    "pressure_temperature",  # good for simple gas tanks
    "temperature_mass",
]

@dataclass
class TankInitialCondition:
    mode: TankInitMode
    pressure_pa: float | None = None
    total_mass_kg: float | None = None
    temperature_k: float | None = None

    def __post_init__(self):
        if self.mode == "pressure_mass":
            if self.pressure_pa is None or self.total_mass_kg is None:
                raise ValueError("pressure_mass requires pressure_pa and total_mass_kg")
        
        elif self.mode == "pressure_temperature":
            if self.pressure_pa is None or self.temperature_k is None:
                raise ValueError("pressure_temperature requires pressure_pa and temperature_k")

        elif self.mode == "temperature_mass":
            if self.temperature_k is None or self.total_mass_kg is None:
                raise ValueError("temperature_mass requires temperature_k and total_mass_kg")

@dataclass
class TankConfig:
    name: str               # eg "ethanol tank", "n2o tank"
    role: SourceRole        # fuel / oxidiser / pressurant
    fluid: Fluid
    tank_volume_m3: float
    phase_model: PhaseModel = 'unknown'
    state: TankState | None = None

    def calculate_internal_energy(
            self,
            total_mass_kg: float,
            temperature_k: float
    ) -> float:
        """
        Find the total internal energy of a two-phase saturated state with the specified total mass and temperature.
        """

        # find the saturation properties at the given temperature
        saturation_properties = self.fluid.get_saturation_properties_from_temp(temperature_k)

        vf = saturation_properties['vf']
        vg = saturation_properties['vg']
        uf = saturation_properties['uf']
        ug = saturation_properties['ug']

        # calculate the mass split at the given temperature

        mass_vapour_kg = (total_mass_kg*vf - self.tank_volume_m3) / (vf - vg)
        mass_liquid_kg = total_mass_kg - mass_vapour_kg

        # calculate the internal energy at the given temperature
        total_internal_energy_j = mass_liquid_kg * uf + mass_vapour_kg * ug


        return total_internal_energy_j
    


    def _initialise_self_pressurised_tank_from_pressure_mass(
            self, 
            initial_condition: TankInitialCondition
    ) -> TankState:
        
        """
        Initialise a saturated self-pressurised tank from total mass and tank pressure.

        Assumes the tank is in liquid-vapour equilibrium and that the specified
        pressure lies between the fluid triple and critical pressures.
        """



        if initial_condition.pressure_pa is None or initial_condition.total_mass_kg is None:
            raise ValueError("pressure_mass requires pressure_pa and total_mass_kg")

        pressure_pa = initial_condition.pressure_pa
        total_mass_kg = initial_condition.total_mass_kg

        # find the saturation temperature given the inital pressure
        saturation_properties = self.fluid.get_saturation_properties_from_pressure(pressure_pa)
        saturation_temp_k = saturation_properties['tsat']

        vf = saturation_properties['vf']    # liquid specific volume (m3 / kg)
        vg = saturation_properties['vg']    # vapour specific volume (m3 / kg)


        mass_vapour_kg = (self.tank_volume_m3 - total_mass_kg*vf)/ (vg - vf)
        mass_liquid_kg = total_mass_kg - mass_vapour_kg

        # validate the phase split
        if mass_vapour_kg < 0.0 or mass_vapour_kg > total_mass_kg:
            raise ValueError(f"Initial condition for {self.name} is not consistent with a two-phase saturated state at {pressure_pa} Pa.")

        liquid_volume_m3 = mass_liquid_kg * vf
        ullage_volume_m3 = self.tank_volume_m3 - liquid_volume_m3   
        
        # validate the volumes
        vapour_volume_m3 = mass_vapour_kg * vg
        volume_error = abs(vapour_volume_m3 - ullage_volume_m3)
        if volume_error > 1e-7:
            raise ValueError(f"Volume calculation failed. Volume error of {volume_error} m^3.")

        # calculate the total internal energy
        liquid_energy_j = mass_liquid_kg * saturation_properties['uf']
        vapour_energy_j = mass_vapour_kg * saturation_properties['ug']
        total_internal_energy_j = liquid_energy_j + vapour_energy_j



        return TankState(
            config=self,
            pressure_pa=pressure_pa,
            temperature_k=saturation_temp_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=mass_liquid_kg,
            vapour_mass_kg=mass_vapour_kg,
            pressurant_gas_mass_kg=None,
            ullage_volume_m3=ullage_volume_m3,
            liquid_volume_m3=liquid_volume_m3
        )




    def initialise_tank_state(
        self,
        initial_condition: TankInitialCondition
    ) -> TankState:

        if self.phase_model == "self_pressurised" and initial_condition.mode == "pressure_mass":
            return self._initialise_self_pressurised_tank_from_pressure_mass(initial_condition)

        raise NotImplementedError(f"Tank initialisation not implemented for phase_model={self.phase_model!r}, mode={initial_condition.mode!r}")


    def state_from_mass_and_energy(
        self,
        total_mass_kg: float,
        total_internal_energy_j: float,
        previous_state: TankState | None = None
    ) -> TankState:

        if previous_state is not None and previous_state.temperature_k is not None:
            guessed_temperature_k = previous_state.temperature_k
        else:
            guessed_temperature_k = self.fluid.get_Ttriple() + 10.0   # arbitrary guess above the triple point

        # find the total internal energy at the guessed energy

        guessed_energy_j = self.calculate_internal_energy(total_mass_kg, guessed_temperature_k)

        # now need to root search to find the temperature which gives the correct internal energy for the given mass, using a bisection method

        T_low_k = self.fluid.get_Ttriple()
        T_high_k = self.fluid.get_Tcrit()

        energy_error_j = guessed_energy_j - total_internal_energy_j
        iteration = 0

        while abs(energy_error_j) > 1e-3:
            if energy_error_j > 0.0:
                T_high_k = guessed_temperature_k
            else:
                T_low_k = guessed_temperature_k
            
            guessed_temperature_k = 0.5 * (T_low_k + T_high_k)
            guessed_energy_j = self.calculate_internal_energy(total_mass_kg, guessed_temperature_k)
            energy_error_j = guessed_energy_j - total_internal_energy_j
            iteration += 1
        
            if iteration > 100:
                raise RuntimeError(f"Root search failed to converge after 100 iterations. Final energy error of {energy_error_j} J.")
        
        # at this point guessed_temperature_k should be the temperature which gives the correct internal energy for the given mass, so we can calculate the rest of the state properties at this temperature

        correct_temperature_k = guessed_temperature_k

        saturation_properties = self.fluid.get_saturation_properties_from_temp(correct_temperature_k)

        vf = saturation_properties['vf']
        vg = saturation_properties['vg']

        mass_vapour_kg = (self.tank_volume_m3 - total_mass_kg*vf)/ (vg - vf)
        mass_liquid_kg = total_mass_kg - mass_vapour_kg

        liquid_volume_m3 = mass_liquid_kg * vf
        vapour_volume_m3 = mass_vapour_kg * vg
        ullage_volume_m3 = self.tank_volume_m3 - mass_liquid_kg * vf


        # validate the volumes
        volume_error = abs(vapour_volume_m3 - ullage_volume_m3)
        if volume_error > 1e-7:
            raise ValueError(f"Volume calculation failed. Volume error of {volume_error} m^3.")


        liquid_energy_j = mass_liquid_kg * saturation_properties['uf']
        vapour_energy_j = mass_vapour_kg * saturation_properties['ug']
        total_internal_energy_j = liquid_energy_j + vapour_energy_j



        return TankState(
            config=self,
            pressure_pa=saturation_properties['psat'],
            temperature_k=correct_temperature_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=mass_liquid_kg,
            vapour_mass_kg=mass_vapour_kg,
            pressurant_gas_mass_kg=None,
            ullage_volume_m3=ullage_volume_m3,
            liquid_volume_m3=liquid_volume_m3
            )


@dataclass
class TankState:
    config: TankConfig

    pressure_pa: float | None = None
    temperature_k: float | None = None
    
    total_mass_kg: float | None = None
    total_internal_energy_j: float | None = None

    liquid_mass_kg: float | None = None
    vapour_mass_kg: float | None = None
    pressurant_gas_mass_kg: float | None = None

    ullage_volume_m3: float | None = None
    liquid_volume_m3: float | None = None




