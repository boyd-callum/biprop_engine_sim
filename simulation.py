
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


# ------------------------------------
# Functions to advance the simulation timestep by one timestep




def gas_advance_timestep(
    tank_config: TankConfig,
    injector_config: InjectorConfig,
    dt_s: float
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

    downstream_pressure = ATMOSPHERE_PRESSURE_PA


    # find current bulk gas density from total mass and tank volume
    density_kg_m3 = tank_config.state.total_mass_kg / tank_config.tank_volume_m3


    # use this calculated density to find the mdot through the injector
    injector_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
        tank_state=tank_config.state,
        downstream_pressure_pa=downstream_pressure
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

    downstream_pressure_pa=ATMOSPHERE_PRESSURE_PA
    # downstream_pressure_pa = tank_pressure_pa * 0.36 # emulate chamber pressure of about 20 bar

    # find the saturation properties at the current tank temperature, to determine the density of the liquid phase in the tank
    saturation_properties = tank_config.fluid.get_saturation_properties_from_temp(tank_temperature_k)
    

    if tank_config.state.liquid_mass_kg is None:
        raise ValueError("Tank state is missing liquid mass.")
    
    liquid_mass_kg = tank_config.state.liquid_mass_kg


    # if any liquid remains at start of time step, estimate mdot via the dyer model

    if liquid_mass_kg > 1e-9:

        # find the saturation properties at the current tank temperature, to determine the density of the liquid phase in the tank
        saturation_properties = tank_config.fluid.get_saturation_properties_from_temp(tank_temperature_k)

        prelim_mdot_kg_s = injector_config.get_dyer_mdot_kg_s(
            tank_state=tank_config.state,
            downstream_pressure_pa=downstream_pressure_pa
        )

        if liquid_mass_kg > prelim_mdot_kg_s * dt_s:
            # entire timestep is liquid discharge through injector
            injector_mdot_kg_s = prelim_mdot_kg_s
            mass_out_kg = injector_mdot_kg_s*dt_s
            energy_out_j = mass_out_kg * saturation_properties["hf"]

        else:
            # dryout occurs in this timestep
            # flow out is still vapour, but the tank state at the start of the timestep is still saturated, so using saturated vapour enthalpy here

            injector_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
                tank_state=tank_config.state,
                downstream_pressure_pa=downstream_pressure_pa
            )
            mass_out_kg = injector_mdot_kg_s*dt_s
            energy_out_j = mass_out_kg*saturation_properties["hg"]

    
    else:
        # tank was already gas-only at the start of the timestep

        injector_mdot_kg_s = injector_config.get_gas_mdot_kg_s(
            tank_state=tank_config.state,
            downstream_pressure_pa=downstream_pressure_pa
        )
        mass_out_kg = injector_mdot_kg_s*dt_s

        # for a real gas-only state we aren't saturated, so need to use actual gas enthalpy
        gas_density_kg_m3 = tank_config.state.total_mass_kg / tank_config.tank_volume_m3
        gas_enthalpy_j_kg = tank_config.fluid.props_si(
            "H", "D", gas_density_kg_m3, "T", tank_temperature_k
        )

        energy_out_j = mass_out_kg * gas_enthalpy_j_kg
   

    # find new total mass and internal energy in the tank
    new_total_mass_kg = tank_config.state.total_mass_kg - mass_out_kg
    new_total_internal_energy_j = tank_config.state.total_internal_energy_j - energy_out_j


    # decide which tank-state solver to use for the new state
    # if liquid still remains, stay on self-pressurised 2phase solver, otherwise switch to gas
    liquid_remains_after_step = (
        tank_config.phase_model == "self_pressurised"
        and tank_config.state.liquid_mass_kg is not None
        and tank_config.state.liquid_mass_kg > mass_out_kg + 1e-9
    )

    

    # update tank state based on the new mass and energy
    new_tank_state = tank_config.state_from_mass_and_energy(
        total_mass_kg=new_total_mass_kg,
        total_internal_energy_j=new_total_internal_energy_j,
        previous_state=tank_config.state,
        phase_override="self_pressurised" if liquid_remains_after_step else "single_phase"
    )

    return new_tank_state, injector_mdot_kg_s
    




#---------------------------
# Overall simulation functions

def gas_simulate(
    case:SimCase,
    record: bool = True
    ) -> SimRecord | None:
    
    """
    Simulate a cold-gas test of a single gas tank with a single injector, using the given case.
    """

    if case.tank_configs is None or len(case.tank_configs) == 0:
        raise ValueError("No tank configs specified in case.")
    if case.injector_configs is None or len(case.injector_configs) == 0:
        raise ValueError("No injector configs specified in case.")


    tank_config = case.tank_configs["n2_tank"]
    injector_config = case.injector_configs["test_injector"]
    inital_condition = case.tank_initial_conditions["n2_tank"]

    # initialise tank state from the given initial condition
    tank_config.state = tank_config.initialise_tank_state(inital_condition)

    dt_s = case.settings.dt_s
    t_final_s = case.settings.t_final_s

    time_s = 0.0
    sim_points = []

    # record the inital state
    sim_points.append(SimPoint(
        time_s=time_s,
        engine=None,
        tanks={"n2_tank" : tank_config.state},
        injectors_mdot={'test_injector' : 0.0}
    ))

    # while loop to advance the simulation in time
    while time_s < t_final_s:

        # advance the simulation by one timestep
        try:
            new_tank_state = gas_advance_timestep(
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
            tanks={'n2_tank': tank_config.state},
            injectors_mdot={'test_injector': injector_mdot_kg_s}
        ))


        # end simulation if tank pressure drops to 110% of atmo pressure
        if tank_config.state.pressure_pa is not None and tank_config.state.pressure_pa <= ATMOSPHERE_PRESSURE_PA*1.5:
            
            print(f"Tank pressure equals 110% of atmospheric at time {time_s:.2f} s, ending simulation.")
            break
        
        # print data for debugging
        if tank_config.state.pressure_pa is not None:
            print(f"time: {time_s:.2f} s,   pressure: {tank_config.state.pressure_pa/1e5:.2f} bar,  total mass: {tank_config.state.total_mass_kg:.2f} kg,    mdot: {injector_mdot_kg_s:.3f} kg")

    if record:
        return SimRecord(points=sim_points)
    else:
        return None




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
        if tank_config.state.pressure_pa is not None and tank_config.state.pressure_pa <= ATMOSPHERE_PRESSURE_PA*1.5:
            
            print(f"Tank pressure equals 110% of atmospheric at time {time_s:.2f} s, ending simulation.")
            break
        
        # print data for debugging
        if tank_config.state.pressure_pa is not None:
            print(f"time: {time_s:.2f} s,   pressure: {tank_config.state.pressure_pa/1e5:.2f} bar,  total mass: {tank_config.state.total_mass_kg:.2f} kg,    massout: {injector_mdot_kg_s*dt_s:.3f} kg")

    if record:
        return SimRecord(points=sim_points)
    else:
        return None
    



def single_tank_simulate(
    case: SimCase,
    record: bool = True
) -> SimRecord | None:
    
    """
    Generic simulator for a single tank blowdown case, with a single injector
    """

    if case.tank_configs is None or len(case.tank_configs) == 0:
        raise ValueError(f"No tank configs included in case {case.name}.")
    if case.injector_configs is None or len(case.injector_configs) == 0:
        raise ValueError(f"No injector configs included in case {case.name}.")
    
    if len(case.tank_configs) > 1:
        raise ValueError(f"single_tank_simulate is only for simulating one tank. {len(case.tank_configs)} tank configs where given.")
    if len(case.injector_configs) > 1:
        raise ValueError(f"single_tank_simulate is only for simulating one injector. {len(case.injector_configs)} injector configs where given.")
    if len(case.tank_initial_conditions) > 1:
        raise ValueError(f"single_tank_simulate is only for simulating one tank. {len(case.injector_configs)} tank initial conditions where given.")


    # get the first entry in both tank and injector congifs. as single tank / injector, should be only entry
    tank_config = list(case.tank_configs.values())[0]
    injector_config = list(case.injector_configs.values())[0]
    initial_condition = list(case.tank_initial_conditions.values())[0]

    # initialise the tank state from the initial condition
    tank_config.state = tank_config.initialise_tank_state(initial_condition)

    dt_s = case.settings.dt_s
    t_final_s = case.settings.t_final_s

    time_s = 0.0
    sim_points = []

    # save the initial state to the record
    sim_points.append(SimPoint(
        time_s = time_s,
        engine = None,
        tanks = {'tank' : tank_config.state},
        injectors_mdot={'injector' : 0.0}
    ))

    # While loop to advance the simulation in time
    while time_s < t_final_s:
        
        # Advance the simulation by one timestep
        try:

            if tank_config.phase_model == "self_pressurised":
                new_tank_state = blowdown_advance_timestep(
                    tank_config=tank_config,
                    injector_config=injector_config,
                    dt_s=dt_s
                )
            if tank_config.phase_model == "single_phase":
                new_tank_state = gas_advance_timestep(
                    tank_config=tank_config,
                    injector_config=injector_config,
                    dt_s=dt_s
                )
            else:
                raise NotImplementedError(f"Phase model {tank_config.phase_model} has not yet been implmented.")

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
            tanks={'tank': tank_config.state},
            injectors_mdot={'injector': injector_mdot_kg_s}
        ))


        # end simulation if tank pressure drops to 110% of atmo pressure
        if tank_config.state.pressure_pa is not None and tank_config.state.pressure_pa <= ATMOSPHERE_PRESSURE_PA*1.5:
            
            print(f"Tank pressure equals 110% of atmospheric at time {time_s:.2f} s, ending simulation.")
            break
        
        # print data for debugging
        if tank_config.state.pressure_pa is not None:
            print(f"time: {time_s:.2f} s,   pressure: {tank_config.state.pressure_pa/1e5:.2f} bar,  total mass: {tank_config.state.total_mass_kg:.2f} kg,    massout: {injector_mdot_kg_s*dt_s:.3f} kg")

    if record:
        return SimRecord(points=sim_points)
    else:
        return None