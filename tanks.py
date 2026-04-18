

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
    injector: InjectorConfig | None = None


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


        return TankState(
            config=self,
            pressure_pa=pressure_pa,
            temperature_k=saturation_temp_k,
            total_mass_kg=total_mass_kg,
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







@dataclass
class TankState:
    config: TankConfig

    pressure_pa: float | None = None
    temperature_k: float | None = None
    
    total_mass_kg: float | None = None 
    liquid_mass_kg: float | None = None
    vapour_mass_kg: float | None = None
    pressurant_gas_mass_kg: float | None = None

    ullage_volume_m3: float | None = None
    liquid_volume_m3: float | None = None

    mass_flow_out_kg_s: float | None = None
    injector_pressure_drop_pa: float | None = None


