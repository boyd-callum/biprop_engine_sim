
from dataclasses import dataclass


from tanks import TankConfig, TankInitialCondition
from fluid import Fluid
from injector import InjectorConfig
from engine import EngineConfig


#------------------------------
# dataclasses

@dataclass
class SimulationSettings:
    dt_s: float = 0.01
    t_final_s: float = 10.0


@dataclass
class SimCase:
    name: str

    # hardware config
    tank_configs: dict[str, TankConfig]
    tank_initital_conditions: dict[str, TankInitialCondition]
    engine_config: EngineConfig

    # sim settings
    settings: SimulationSettings

#------------------------------
# defining fluids

NITROUS_OXIDE = Fluid(
    name="nitrous_oxide",
    coolprop_name="NitrousOxide",
    backend="HEOS",
    cea_name="N2O",
)

ETHANOL = Fluid(
    name="ethanol",
    coolprop_name="Ethanol",
    backend="HEOS",
    cea_name="ETHANOL",
)

NITROGEN = Fluid(
    name="nitrogen",
    coolprop_name="Nitrogen",
    backend="HEOS",
    cea_name=None,   # not used in combustion
)


#------------------------------
# N2O Blowdown case

