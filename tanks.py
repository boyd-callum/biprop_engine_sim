

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
    temperature_k: float | None = None
    total_mass_kg: float | None = None

@dataclass
class TankConfig:
    name: str               # eg "ethanol tank", "n2o tank"
    role: SourceRole        # fuel / oxidiser / pressurant
    fluid: Fluid
    tank_volume_m3: float
    phase_model: PhaseModel = 'unknown'
    injector: InjectorConfig | None = None


@dataclass
class TankState:
    config: TankConfig

    pressure_pa: float | None = None
    temperature_k: float | None = None
    
    total_mass_kg: float | None = None      # derived from phase masses
    liquid_mass_kg: float | None = None
    vapour_mass_kg: float | None = None
    pressurant_gas_mass_kg: float | None = None

    ullage_volume_m3: float | None = None
    liquid_volume_m3: float | None = None

    mass_flow_out_kg_s: float | None = None
    injector_pressure_drop_pa: float | None = None



   