
from dataclasses import dataclass, field


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
    tank_initial_conditions: dict[str, TankInitialCondition]

    injector_configs: dict[str, InjectorConfig]

    engine_config: EngineConfig | None = None

    # sim settings
    settings: SimulationSettings = field(default_factory=SimulationSettings)

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

test_injector = InjectorConfig(
    cd=0.75,
    area_m2=5e-5
)

n2o_tank = TankConfig(
    name="Nitrous Oxide Tank",
    role="oxidiser",
    fluid=NITROUS_OXIDE,
    tank_volume_m3=0.01,    # 10L
    phase_model="self_pressurised"
)

n2o_tank_initial = TankInitialCondition(
    mode="pressure_mass",
    pressure_pa=6e+6,    # 60 bar
    total_mass_kg=6.0
)  

#------------------------------
# N2O Blowdown case

n2o_blowdown_case = SimCase(
    name="n2o_blowdown_test",
    tank_configs={"n2o_tank": n2o_tank},
    tank_initial_conditions={"n2o_tank": n2o_tank_initial},
    injector_configs={"test_injector":test_injector},
    engine_config=None,
    settings=SimulationSettings(0.01, 20)
)