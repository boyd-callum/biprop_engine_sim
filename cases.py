
from dataclasses import dataclass, field


from tanks import TankConfig, TankInitialCondition
from fluid import Fluid
from injector import InjectorConfig
from engine import EngineConfig
from regulator import RegulatorConfig

from constants import ATMOSPHERE_PRESSURE_PA

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
    name="Nitrogen Tank",
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
    name="Ethanol Tank",
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
    pressure_pa=60*1e5 # 60 bar
)

ethanol_regulator = RegulatorConfig(
    name="Ethanol Reg",
    set_pressure_pa=60*1e5 # 60 bar
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




if __name__ == "__main__":

    case = n2_blowdown_case

    tank_config = case.tank_configs["n2_tank"]
    injector_config = case.injector_configs["test_injector"]
    inital_condition = case.tank_initial_conditions["n2_tank"]

    # initialise tank state from the given initial condition
    tank_config.state = tank_config.initialise_tank_state(inital_condition)


    print(tank_config.state)

    injector_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
        tank_state=tank_config.state,
        downstream_pressure_pa=ATMOSPHERE_PRESSURE_PA
    )

    print(injector_mdot_kg_s)

    total_mass_kg = tank_config.state.total_mass_kg
    total_internal_energy_j = tank_config.state.total_internal_energy_j

    if total_internal_energy_j is None or total_mass_kg is None:
        raise ValueError("mass or energy are None")

    # use tank state solver instead of initialise solver
    tank_config.state = tank_config.state_from_mass_and_energy(
        total_mass_kg=total_mass_kg,
        total_internal_energy_j=total_internal_energy_j,
        phase_override="single_phase"
    )

    print(tank_config.state)
