from __future__ import annotations

from dataclasses import dataclass, field


from tanks import TankConfig, TankInitialCondition
from fluid import Fluid
from injector import InjectorConfig
from engine import EngineConfig, EngineGeometry
from regulator import RegulatorConfig
from openrocket_exporter import RSEComponent

from constants import ATMOSPHERE_PRESSURE_PA, PI

#------------------------------
# dataclasses

@dataclass
class SimulationSettings:
    dt_s: float = 0.01
    t_final_s: float = 10.0
    print_updates: bool = True
    print_steps: int = 10


@dataclass
class SimCase:
    name: str

    # hardware config
    tank_configs: dict[str, TankConfig]
    tank_initial_conditions: dict[str, TankInitialCondition]

    injector_configs: dict[str, InjectorConfig]

    regulator_configs: dict[str, RegulatorConfig] | None = None

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


# ----------------------
# Defining configs


test_injector = InjectorConfig(
    cd=0.65,
    area_m2=6.234e-5, # 8 holes of area 3.15mm
    k=2 # approximation
)

n2_test_injector = InjectorConfig(
    cd=0.65,
    area_m2=1e-5,
    k=2 # approximation
)

n2o_tank = TankConfig(
    name="Nitrous Oxide Tank",
    role="oxidiser",
    fluid=NITROUS_OXIDE,
    tank_volume_m3=0.0088,    # 8.8L
    phase_model="self_pressurised"
)

n2o_tank_initial = TankInitialCondition(
    mode="pressure_mass",
    pressure_pa=60e+5,    # 60 bar
    total_mass_kg=5.85
)  

n2_tank = TankConfig(
    name="n2_tank",
    role="pressurant",
    fluid=NITROGEN,
    tank_volume_m3=0.001,    # 1L
    phase_model="single_phase"
)

n2_tank_initial = TankInitialCondition(
    mode="pressure_temperature",
    pressure_pa=3e+7,   # 300 bar
    temperature_k=293.15 # 20C
)


ethanol_tank = TankConfig(
    name="ethanol_tank",
    role="fuel",
    fluid=ETHANOL,
    tank_volume_m3=0.005, # 5 L
    phase_model="pressurised_liquid",
    pressurant_fluid=NITROGEN
)

ethanol_tank_initial = TankInitialCondition(
    mode="pressure_temperature_mass",
    total_mass_kg=3, # actually liquid mass
    temperature_k=293.15, # 20c
    pressure_pa=50*1e5 # 60 bar
)

ethanol_regulator = RegulatorConfig(
    name="ethanol_regulator",
    set_pressure_pa=50*1e5, # 50 bar
    role="fuel"
)


throat_dia_mm = 25
throat_area_mm2 = PI * (throat_dia_mm/2)**2
throat_area_m2 = throat_area_mm2 * 1e-6

engine_geometry = EngineGeometry(
    nozzle_throat_area_m2 = throat_area_m2,
    expansion_ratio = 6.8
)

engine = EngineConfig(
    geometry = engine_geometry
)


n2o_injector = InjectorConfig(
    cd=0.65,
    area_m2=3.0e-5,
    role="oxidiser",
)

ethanol_injector = InjectorConfig(
    cd=0.65,
    area_m2=6.5e-6,
    role="fuel",
)


#------------------------------
# Squirtle sizing

squirtle_n2_tank = TankConfig(
    name="squirtle_nitrogen_tank",
    role="pressurant",
    fluid=NITROGEN,
    tank_volume_m3=0.00075,    # 0.75L
    phase_model="single_phase"
)
squirtle_n2_tank_initial = TankInitialCondition(
    mode="pressure_temperature",
    pressure_pa=3.5e+7,   # 350 bar
    temperature_k=293.15 # 20C
)


squirtle_ethanol_tank = TankConfig(
    name="squirtle_ethanol_tank",
    role="fuel",
    fluid=ETHANOL,
    tank_volume_m3=0.0022, # 2.2 L
    phase_model="pressurised_liquid",
    pressurant_fluid=NITROGEN
)
squirtle_ethanol_tank_initial = TankInitialCondition(
    mode="pressure_temperature_mass",
    total_mass_kg=1.5, # actually liquid mass
    temperature_k=293.15, # 20c
    pressure_pa=50*1e5 # 50 bar
)
squirtle_ethanol_regulator = RegulatorConfig(
    name="squirtle_ethanol_regulator",
    set_pressure_pa=50*1e5, # 50 bar
    role="fuel"
)


squirtle_n2o_tank = TankConfig(
    name="squirtle_n2o_tank",
    role="oxidiser",
    fluid=NITROUS_OXIDE,
    tank_volume_m3=0.01,    # 10L
    phase_model="self_pressurised"
)
squirtle_n2o_tank_initial = TankInitialCondition(
    mode="pressure_mass",
    pressure_pa=60e+5,    # 60 bar
    total_mass_kg=5
)


squirtle_n2o_injector = InjectorConfig(
    cd=0.75,
    area_m2=1.77e-5, # 15.7mm^2
    role="oxidiser"
)
squirtle_ethanol_injector = InjectorConfig(
    cd=0.75,
    area_m2=3.17e-6, # 3mm^2
    role="fuel"
)



squirtle_engine_geometry = EngineGeometry(
    nozzle_throat_area_m2 = 4.15e-4, # throat diameter 23mm
    expansion_ratio = 4.3
)
squirtle_engine = EngineConfig(
    geometry = squirtle_engine_geometry,
    cstar_efficiency=0.90,
    cf_efficiency=0.95
)



#------------------------------
# Cases

n2o_blowdown_case = SimCase(
    name="n2o_blowdown_test",
    tank_configs={"n2o_tank": n2o_tank},
    tank_initial_conditions={"n2o_tank": n2o_tank_initial},
    injector_configs={"test_injector":test_injector},
    engine_config=None,
    settings=SimulationSettings(0.01, 20)
)



n2_blowdown_case = SimCase(
    name="n2_blowdown_test",
    tank_configs={"n2_tank": n2_tank},
    tank_initial_conditions={"n2_tank": n2_tank_initial},
    injector_configs={"test_injector":n2_test_injector},
    engine_config=None,
    settings=SimulationSettings(0.01, 20)
)

ethanol_case = SimCase(
    name="ethanol_test",
    tank_configs={
            "ethanol_tank": ethanol_tank,
            "n2_tank": n2_tank
        },
    tank_initial_conditions={
            "ethanol_tank": ethanol_tank_initial,
            "n2_tank": n2_tank_initial
        },
    injector_configs={"test_injector":test_injector},
    regulator_configs={"ethanol_regulator": ethanol_regulator},
    engine_config=None,
    settings=SimulationSettings(0.01, 20)
)




full_biprop_case = SimCase(
    name = "Full Biprop Case",

    tank_configs = {
        "n2o_tank" : n2o_tank,
        "ethanol_tank" : ethanol_tank,
        "n2_tank" : n2_tank
        },

    tank_initial_conditions={
        "n2o_tank" : n2o_tank_initial,
        "ethanol_tank": ethanol_tank_initial,
        "n2_tank": n2_tank_initial
        },

    injector_configs={
        "n2o_injector" : n2o_injector,
        "ethanol_injector" : ethanol_injector
        },

    regulator_configs={
        "ethanol_regulator" : ethanol_regulator
        },
    
    engine_config=engine,

    settings=SimulationSettings(dt_s=0.02, print_steps=10)

)


squirtle_case = SimCase(
    name = "Squirtle Biprop Case",

    tank_configs = {
        "n2o_tank" : squirtle_n2o_tank,
        "ethanol_tank" : squirtle_ethanol_tank,
        "n2_tank" : squirtle_n2_tank
        },

    tank_initial_conditions={
        "n2o_tank" : squirtle_n2o_tank_initial,
        "ethanol_tank": squirtle_ethanol_tank_initial,
        "n2_tank": squirtle_n2_tank_initial
        },

    injector_configs={
        "n2o_injector" : squirtle_n2o_injector,
        "ethanol_injector" : squirtle_ethanol_injector
        },

    regulator_configs={
        "ethanol_regulator" : squirtle_ethanol_regulator
        },
    
    engine_config=squirtle_engine,

    settings=SimulationSettings(dt_s=0.02, t_final_s=15, print_steps=50, print_updates=True)

)



#-------------------------------
# OpenRocket exports

