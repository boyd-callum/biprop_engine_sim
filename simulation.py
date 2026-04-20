
from dataclasses import dataclass
from typing import Literal


from injector import InjectorConfig
from tanks import TankState, TankConfig
from engine import EngineState
from cases import SimCase
from constants import *

# -------------------------------
# Dataclasses for the simulation running

 


@dataclass
class SimPoint:
    """A point in the simulation, with a timestamp and states for the engine and fluid sources."""
    time_s: float

    engine: EngineState | None = None
    tanks: dict[str, TankState] | None = None
    injectors_mdot: dict[str, float] | None = None
    
@dataclass
class SimRecord:
    """A record of the simulation, containing a list of SimPoints."""
    points: list[SimPoint]


def blowdown_advance_timestep(
     tank_config: TankConfig,
     injector_config: InjectorConfig,
     dt_s: float
    ) -> tuple[TankState, float]:

    if tank_config.state is None:
        raise ValueError("Tank state is not initialised.")
    if tank_config.state.pressure_pa is None or tank_config.state.temperature_k is None:
        raise ValueError("Tank state is missing pressure or temperature.")
    if tank_config.state.total_mass_kg is None or tank_config.state.total_internal_energy_j is None:
        raise ValueError("Tank state is missing total mass or total internal energy.")

    # find current tank pressure and temperature

    tank_pressure_pa = tank_config.state.pressure_pa
    tank_temperature_k = tank_config.state.temperature_k

    # find the saturation properties at the current tank temperature, to determine the density of the liquid phase in the tank

    saturation_properties = tank_config.fluid.get_saturation_properties_from_temp(tank_temperature_k)

    liquid_density_kg_m3 = 1.0/saturation_properties['vf']

    # find the mass flow rate through the injector at the current tank pressure and temperature

    injector_mdot_kg_s = injector_config.get_liquid_mdot_kg_s(
        upstream_pressure_pa=tank_pressure_pa,
        downstream_pressure_pa=AMTOSPHERE_PRESSURE_PA,  # assume atmospheric pressure downstream
        liquid_density_kg_m3=liquid_density_kg_m3
    )

    # find the mass of fluid leaving the tank during this timestep
    mass_out_kg = injector_mdot_kg_s * dt_s

    # find energy leaving the tank during this timestep, assuming the fluid leaving is saturated liquid at the tank temperature
    energy_out_j = mass_out_kg * saturation_properties['uf']

    # find new total mass and internal energy in the tank
    new_total_mass_kg = tank_config.state.total_mass_kg - mass_out_kg
    new_total_internal_energy_j = tank_config.state.total_internal_energy_j - energy_out_j

    # update tank state based on the new mass and energy
    new_tank_state = tank_config.state_from_mass_and_energy(
        total_mass_kg=new_total_mass_kg,
        total_internal_energy_j=new_total_internal_energy_j
    )

    return new_tank_state, injector_mdot_kg_s
    


def blowdown_simulate(
        case: SimCase,
        record: bool = True
) -> SimRecord | None:  
    
    """
    Simulate a cold-gas blowdown of a single self-pressurised tank with a single injector, using the specified case.
    """

    if case.tank_configs is None or len(case.tank_configs) == 0:
        raise ValueError("No tank configs specified in case.")
    if case.injector_configs is None or len(case.injector_configs) == 0:
        raise ValueError("No injector configs specified in case.")
    
    tank_config = case.tank_configs['n2o_tank']
    injector_config = case.injector_configs['test_injector']
    initial_condition = case.tank_initial_conditions['n2o_tank']

    # initialise the tank state from the initial condition
    tank_config.state = tank_config.initialise_tank_state(initial_condition)

    dt_s = case.settings.dt_s
    t_final_s = case.settings.t_final_s

    time_s = 0.0
    sim_points = []

    # record the initial state
    sim_points.append(SimPoint(
        time_s=time_s,
        engine=None,
        tanks={'n2o_tank': tank_config.state},
        injectors_mdot={'test_injector': 0.0}
    ))

    # While loop to advance the simulation in time
    while time_s < t_final_s:
        
        # Advance the simulation by one timestep
        try:
            new_tank_state = blowdown_advance_timestep(
                tank_config=tank_config,
                injector_config=injector_config,
                dt_s=dt_s
            )
        except RuntimeError as e:
            print(f"Error advancing timestep at time {time_s:.2f} s: {e}")
            print(f"final state at time of error: pressure={tank_config.state.pressure_pa:.2f} Pa, temperature={tank_config.state.temperature_k:.2f} K, total_mass={tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={tank_config.state.total_internal_energy_j:.2f} J")
            break

        # Update the tank state
        tank_config.state = new_tank_state[0]

        # Update the injector mass flow rate
        injector_mdot_kg_s = new_tank_state[1]

        # Update the time
        time_s += dt_s

        # Record the state at this timestep
        sim_points.append(SimPoint(
            time_s=time_s,
            engine=None,
            tanks={'n2o_tank': tank_config.state},
            injectors_mdot={'test_injector': injector_mdot_kg_s}
        ))


        # end simulation if tank pressure drops to 110% of atmo pressure
        if tank_config.state.pressure_pa is not None and tank_config.state.pressure_pa <= AMTOSPHERE_PRESSURE_PA*1.1:
            
            print(f"Tank pressure equals 110% of atmospheric at time {time_s:.2f} s, ending simulation.")
            break
        
        # print data for debugging
        if tank_config.state.pressure_pa is not None:
            print(f"time: {time_s:.2f} s,   pressure: {tank_config.state.pressure_pa/1e5:.2f} bar")

    if record:
        return SimRecord(points=sim_points)
    else:
        return None
    