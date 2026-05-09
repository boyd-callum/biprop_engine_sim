
from dataclasses import dataclass
from typing import Literal


from injector import InjectorConfig
from tanks import TankState, TankConfig
from engine import EngineState
from cases import SimCase
from regulator import RegulatorConfig
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
    regulators_mdot: dict[str, float] | None = None
    
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

    if liquid_mass_kg > DRYOUT_TOLERANCE_KG:

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
        and tank_config.state.liquid_mass_kg > mass_out_kg + DRYOUT_TOLERANCE_KG
    )

    

    # update tank state based on the new mass and energy
    new_tank_state = tank_config.state_from_mass_and_energy(
        total_mass_kg=new_total_mass_kg,
        total_internal_energy_j=new_total_internal_energy_j,
        previous_state=tank_config.state,
        phase_override="self_pressurised" if liquid_remains_after_step else "single_phase"
    )

    return new_tank_state, injector_mdot_kg_s
    



def pressurised_liquid_advance_timestep(
        liquid_tank_config: TankConfig,
        injector_config: InjectorConfig,
        regulator_config: RegulatorConfig,
        pressurant_tank_config: TankConfig,
        dt_s: float,
        downstream_pressure_pa: float = ATMOSPHERE_PRESSURE_PA
) -> tuple[ TankState, TankState, float, float, bool]:
    

    liquid_state = liquid_tank_config.state
    if liquid_state is None:
        raise ValueError("Liquid tank state is not initialised.")
    if liquid_state.pressure_pa is None or liquid_state.temperature_k is None:
        raise ValueError("Liquid tank state is missing pressure or temperature.")
    if liquid_state.total_mass_kg is None:
        raise ValueError("Liquid tank state is missing total mass.")
    if liquid_state.total_internal_energy_j is None:
        raise ValueError("Liquid tank state is missing total internal energy.")
    if liquid_state.pressurant_gas_mass_kg is None:
        raise ValueError("Liquid tank state is missing pressurant gas mass.")
    if liquid_state.liquid_mass_kg is None:
        raise ValueError("Liquid tank state is missing liquid mass.")
    
    pressurant_state = pressurant_tank_config.state
    if pressurant_state is None:
        raise ValueError("Pressurant tank state is not initialised.")
    if pressurant_state.pressure_pa is None or pressurant_state.temperature_k is None:
        raise ValueError("Pressurant tank state is missing pressure or temperature.")
    if pressurant_state.total_mass_kg is None:
        raise ValueError("Pressurant tank state is missing total mass.")
    if pressurant_state.total_internal_energy_j is None:
        raise ValueError("Pressurant tank state is missing total internal energy.")
    
    if liquid_tank_config.pressurant_fluid is None:
        raise ValueError("Liquid tank config is missing a defined pressurant fluid.")
    
    #-------------------------
    # 1. remove liquid from the liquid tank via the injector

    liquid_pressure_pa = liquid_state.pressure_pa
    liquid_temperature_k = liquid_state.temperature_k
    liquid_mass_kg = liquid_state.liquid_mass_kg
    liquid_pressurant_mass_kg = liquid_state.pressurant_gas_mass_kg
    liquid_tank_total_internal_energy_j = liquid_state.total_internal_energy_j

    liquid_density_kg_m3 = liquid_tank_config.fluid.get_fluid_density_from_pressure_temperature(
        P=liquid_pressure_pa,
        T=liquid_temperature_k
    )

    delta_p_pa = liquid_pressure_pa - downstream_pressure_pa
    if delta_p_pa <= 0.0 or liquid_mass_kg <= 0.0:
        # no flow through the injector
        liquid_mdot_kg_s = 0.0

    else:
        liquid_mdot_kg_s = injector_config.get_liquid_mdot_kg_s(
            tank_state=liquid_state,
            downstream_pressure_pa=downstream_pressure_pa
        )
    
    liquid_mass_out_kg = liquid_mdot_kg_s * dt_s

    # clamp so we never remove more liquid mass than is actually in the tank
    if liquid_mass_out_kg > liquid_mass_kg:
        liquid_mass_out_kg = liquid_mass_kg
        liquid_mdot_kg_s = liquid_mass_out_kg / dt_s
    
    # for energy out, we can use the enthalpy of the liquid at current tank pressure and temperature, as the flow is liquid at the injector inlet
    liquid_enthalpy_j_kg = liquid_tank_config.fluid.get_fluid_enthalpy_from_pressure_temperature(
        P=liquid_pressure_pa,
        T=liquid_temperature_k
    )
    liquid_energy_out_j = liquid_mass_out_kg * liquid_enthalpy_j_kg

    new_liquid_mass_kg = liquid_mass_kg - liquid_mass_out_kg

    liquid_dryout = new_liquid_mass_kg <= DRYOUT_TOLERANCE_KG
    if liquid_dryout:
        new_liquid_mass_kg = 0.0

    
    #-------------------------
    # 2. add persurant through the regulator

    pressurant_mdot_kg_s = regulator_config.get_regulator_mdot_kg_s(
        liquid_tank_state=liquid_state,
        pressurant_tank_state=pressurant_state,
        new_liquid_mass_kg=new_liquid_mass_kg,
        dt_s=dt_s
    )

    pressurant_mass_in_kg = pressurant_mdot_kg_s * dt_s

    # clamp againt the available mass in the pressurant tank
    if pressurant_mass_in_kg > pressurant_state.total_mass_kg:
        pressurant_mass_in_kg = pressurant_state.total_mass_kg
        pressurant_mdot_kg_s = pressurant_mass_in_kg / dt_s
    
    pressurant_tank_density_kg_m3 = pressurant_state.total_mass_kg / pressurant_tank_config.tank_volume_m3

    pressurant_inlet_enthalpy_j_kg = liquid_tank_config.pressurant_fluid.get_fluid_enthalpy_from_density_temperature(
        D=pressurant_tank_density_kg_m3,
        T=pressurant_state.temperature_k
    )

    pressurant_energy_in_j = pressurant_mass_in_kg * pressurant_inlet_enthalpy_j_kg

    new_liquid_pressurant_mass_kg = liquid_pressurant_mass_kg + pressurant_mass_in_kg

    #-------------------------
    # 3. update total internal energies

    new_liquid_tank_total_internal_energy_j = liquid_tank_total_internal_energy_j - liquid_energy_out_j + pressurant_energy_in_j

    new_pressurant_tank_total_mass = pressurant_state.total_mass_kg - pressurant_mass_in_kg
    new_pressurant_tank_total_internal_energy_j = pressurant_state.total_internal_energy_j - pressurant_energy_in_j


    if new_pressurant_tank_total_mass <= 0.0:
        raise ValueError("Pressurant tank mass went negative, this should have been clamped against.")
    

    # -------------------------
    # 4. reconstruct tank states

    new_liquid_tank_state = liquid_tank_config.state_from_mass_and_energy(
        liquid_mass_kg=new_liquid_mass_kg,
        pressurant_gas_mass_kg=new_liquid_pressurant_mass_kg,
        total_internal_energy_j=new_liquid_tank_total_internal_energy_j,
        previous_state=liquid_state,
        phase_override="pressurised_liquid",
    )

    new_pressurant_tank_state = pressurant_tank_config.state_from_mass_and_energy(
        total_mass_kg=new_pressurant_tank_total_mass,
        total_internal_energy_j=new_pressurant_tank_total_internal_energy_j,
        previous_state=pressurant_state,
        phase_override="single_phase"
    )

    return (
        new_liquid_tank_state,
        new_pressurant_tank_state,
        liquid_mdot_kg_s,
        pressurant_mdot_kg_s,
        liquid_dryout
    )



#---------------------------
# Overall simulation functions

# def gas_simulate(
#     case:SimCase,
#     record: bool = True
#     ) -> SimRecord | None:
    
#     """
#     Simulate a cold-gas test of a single gas tank with a single injector, using the given case.
#     """

#     if case.tank_configs is None or len(case.tank_configs) == 0:
#         raise ValueError("No tank configs specified in case.")
#     if case.injector_configs is None or len(case.injector_configs) == 0:
#         raise ValueError("No injector configs specified in case.")


#     tank_config = case.tank_configs["n2_tank"]
#     injector_config = case.injector_configs["test_injector"]
#     inital_condition = case.tank_initial_conditions["n2_tank"]

#     # initialise tank state from the given initial condition
#     tank_config.state = tank_config.initialise_tank_state(inital_condition)

#     dt_s = case.settings.dt_s
#     t_final_s = case.settings.t_final_s

#     time_s = 0.0
#     sim_points = []

#     # record the inital state
#     sim_points.append(SimPoint(
#         time_s=time_s,
#         engine=None,
#         tanks={"n2_tank" : tank_config.state},
#         injectors_mdot={'test_injector' : 0.0}
#     ))

#     # while loop to advance the simulation in time
#     while time_s < t_final_s:

#         # advance the simulation by one timestep
#         try:
#             new_tank_state = gas_advance_timestep(
#                 tank_config=tank_config,
#                 injector_config=injector_config,
#                 dt_s=dt_s
#             )
#         except RuntimeError as e:
#             print(f"Error advancing timestep at time {time_s:.2f} s: {e}")
#             print(f"final state at time of error: pressure={tank_config.state.pressure_pa:.2f} Pa, temperature={tank_config.state.temperature_k:.2f} K, total_mass={tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={tank_config.state.total_internal_energy_j:.2f} J")
#             break
        
        
#         # Update the tank state
#         tank_config.state = new_tank_state[0]

#         # Update the injector mass flow rate
#         injector_mdot_kg_s = new_tank_state[1]

#         # Update the time
#         time_s += dt_s

#         # Record the state at this timestep
#         sim_points.append(SimPoint(
#             time_s=time_s,
#             engine=None,
#             tanks={'n2_tank': tank_config.state},
#             injectors_mdot={'test_injector': injector_mdot_kg_s}
#         ))


#         # end simulation if tank pressure drops to 110% of atmo pressure
#         if tank_config.state.pressure_pa is not None and tank_config.state.pressure_pa <= ATMOSPHERE_PRESSURE_PA*1.5:
            
#             print(f"Tank pressure equals 110% of atmospheric at time {time_s:.2f} s, ending simulation.")
#             break
        
#         # print data for debugging
#         if tank_config.state.pressure_pa is not None:
#             print(f"time: {time_s:.2f} s,   pressure: {tank_config.state.pressure_pa/1e5:.2f} bar,  total mass: {tank_config.state.total_mass_kg:.2f} kg,    mdot: {injector_mdot_kg_s:.3f} kg")

#     if record:
#         return SimRecord(points=sim_points)
#     else:
#         return None




# def blowdown_simulate(
#     case: SimCase,
#     record: bool = True
#     ) -> SimRecord | None:  
    
#     """
#     Simulate a cold-gas blowdown of a single self-pressurised tank with a single injector, using the specified case.
#     """

#     if case.tank_configs is None or len(case.tank_configs) == 0:
#         raise ValueError("No tank configs specified in case.")
#     if case.injector_configs is None or len(case.injector_configs) == 0:
#         raise ValueError("No injector configs specified in case.")
    
#     tank_config = case.tank_configs['n2o_tank']
#     injector_config = case.injector_configs['test_injector']
#     initial_condition = case.tank_initial_conditions['n2o_tank']

#     # initialise the tank state from the initial condition
#     tank_config.state = tank_config.initialise_tank_state(initial_condition)

#     dt_s = case.settings.dt_s
#     t_final_s = case.settings.t_final_s

#     time_s = 0.0
#     sim_points = []

#     # record the initial state
#     sim_points.append(SimPoint(
#         time_s=time_s,
#         engine=None,
#         tanks={'n2o_tank': tank_config.state},
#         injectors_mdot={'test_injector': 0.0}
#     ))

#     # While loop to advance the simulation in time
#     while time_s < t_final_s:
        
#         # Advance the simulation by one timestep
#         try:
#             new_tank_state = blowdown_advance_timestep(
#                 tank_config=tank_config,
#                 injector_config=injector_config,
#                 dt_s=dt_s
#             )
#         except RuntimeError as e:
#             print(f"Error advancing timestep at time {time_s:.2f} s: {e}")
#             print(f"final state at time of error: pressure={tank_config.state.pressure_pa:.2f} Pa, temperature={tank_config.state.temperature_k:.2f} K, total_mass={tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={tank_config.state.total_internal_energy_j:.2f} J")
#             break

#         # Update the tank state
#         tank_config.state = new_tank_state[0]

#         # Update the injector mass flow rate
#         injector_mdot_kg_s = new_tank_state[1]

#         # Update the time
#         time_s += dt_s

#         # Record the state at this timestep
#         sim_points.append(SimPoint(
#             time_s=time_s,
#             engine=None,
#             tanks={'n2o_tank': tank_config.state},
#             injectors_mdot={'test_injector': injector_mdot_kg_s}
#         ))


#         # end simulation if tank pressure drops to 110% of atmo pressure
#         if tank_config.state.pressure_pa is not None and tank_config.state.pressure_pa <= ATMOSPHERE_PRESSURE_PA*1.5:
            
#             print(f"Tank pressure equals 110% of atmospheric at time {time_s:.2f} s, ending simulation.")
#             break
        
#         # print data for debugging
#         if tank_config.state.pressure_pa is not None:
#             print(f"time: {time_s:.2f} s,   pressure: {tank_config.state.pressure_pa/1e5:.2f} bar,  total mass: {tank_config.state.total_mass_kg:.2f} kg,    massout: {injector_mdot_kg_s*dt_s:.3f} kg")

#     if record:
#         return SimRecord(points=sim_points)
#     else:
#         return None
    



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
            elif tank_config.phase_model == "single_phase":
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
    




def pressurised_liquid_simulate(
        case: SimCase,
        record: bool = True
) -> SimRecord | None:
    """
    Simulator for a pressurised liquid tank with a single injector, using the given case.
    """
    
    if case.tank_configs is None or len(case.tank_configs) != 2:
        raise ValueError(f"Pressurised liquid simulate requires exactly 2 tank configs, one for the liquid tank and one for the pressurant tank. {len(case.tank_configs) if case.tank_configs is not None else 0} were given.")
    if case.injector_configs is None or len(case.injector_configs) != 1:
        raise ValueError(f"Pressurised liquid simulate requires exactly 1 injector config. {len(case.injector_configs) if case.injector_configs is not None else 0} were given.")
    if case.regulator_configs is None or len(case.regulator_configs) != 1:
        raise ValueError(f"Pressurised liquid simulate requires exactly 1 regulator config. {len(case.regulator_configs) if case.regulator_configs is not None else 0} were given.")
    if case.tank_initial_conditions is None or len(case.tank_initial_conditions) != 2:
        raise ValueError(f"Pressurised liquid simulate requires exactly 2 tank initial conditions, one for the liquid tank and one for the pressurant tank. {len(case.tank_initial_conditions) if case.tank_initial_conditions is not None else 0} were given.")
    

    # determine which tank config is the liquid and which is the pressurant based on the "role" field in the tank config
    liquid_tank_config = None
    pressurant_tank_config = None
    for tank_name, tank_config in case.tank_configs.items():
        if tank_config.role == "fuel" or tank_config.role == "oxidiser":
            liquid_tank_config = tank_config
        elif tank_config.role == "pressurant":
            pressurant_tank_config = tank_config
        else:
            raise ValueError(f"Tank config {tank_name} has invalid role {tank_config.role}. Role must be either 'liquid' or 'pressurant'.")
        
    
    if liquid_tank_config is None:
        raise ValueError("No liquid tank config found in case.")
    if pressurant_tank_config is None:
        raise ValueError("No pressurant tank config found in case.")
    
    # determine which initial condition corresponds to which tank
    liquid_tank_initial_condition = None
    pressurant_tank_initial_condition = None
    for tank_name, initial_condition in case.tank_initial_conditions.items():
        # print(f"tank name: {tank_name}, liquid tank name: {liquid_tank_config.name}, pressurant tank name: {pressurant_tank_config.name}")
        if tank_name == liquid_tank_config.name:
            liquid_tank_initial_condition = initial_condition
        elif tank_name == pressurant_tank_config.name:
            pressurant_tank_initial_condition = initial_condition
        else:
            raise ValueError(f"Tank initial condition {tank_name} does not match any tank config names.")
    
    if liquid_tank_initial_condition is None:
        raise ValueError("No initial condition found for liquid tank.")
    if pressurant_tank_initial_condition is None:
        raise ValueError("No initial condition found for pressurant tank.")


    injector_config = list(case.injector_configs.values())[0]
    regulator_config = list(case.regulator_configs.values())[0]


    # initialise the tank states from the initial conditions
    pressurant_tank_config.state = pressurant_tank_config.initialise_tank_state(pressurant_tank_initial_condition)
    liquid_tank_config.state = liquid_tank_config.initialise_tank_state(liquid_tank_initial_condition)

    # sim settings
    dt_s = case.settings.dt_s
    t_final_s = case.settings.t_final_s

    time_s = 0.0
    sim_points = []

    # record the initial state
    sim_points.append(SimPoint(
        time_s = time_s,
        engine = None,
        tanks = {
            'liquid_tank' : liquid_tank_config.state,
            'pressurant_tank' : pressurant_tank_config.state
        },
        injectors_mdot={'injector' : 0.0}
    ))


    # main while loop to advance the simulation in time
    while time_s < t_final_s:

        # advance the simulation by one timestep
        try:

            (
                new_liquid_tank_state,
                new_pressurant_tank_state,
                liquid_mdot_kg_s,
                pressurant_mdot_kg_s,
                liquid_dryout
            ) = pressurised_liquid_advance_timestep(
                liquid_tank_config=liquid_tank_config,
                injector_config=injector_config,
                regulator_config=regulator_config,
                pressurant_tank_config=pressurant_tank_config,
                dt_s=dt_s
            )
        
        except RuntimeError as e:
            print(f"Error advancing timestep at time {time_s:.2f} s: {e}")
            print(f"final liquid tank state at time of error: pressure={liquid_tank_config.state.pressure_pa:.2f} Pa, temperature={liquid_tank_config.state.temperature_k:.2f} K, total_mass={liquid_tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={liquid_tank_config.state.total_internal_energy_j:.2f} J, liquid_mass={liquid_tank_config.state.liquid_mass_kg:.2f} kg, pressurant_mass={liquid_tank_config.state.pressurant_gas_mass_kg:.2f} kg")
            print(f"final pressurant tank state at time of error: pressure={pressurant_tank_config.state.pressure_pa:.2f} Pa, temperature={pressurant_tank_config.state.temperature_k:.2f} K, total_mass={pressurant_tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={pressurant_tank_config.state.total_internal_energy_j:.2f} J")
            break


        
        # update the tank states
        liquid_tank_config.state = new_liquid_tank_state
        pressurant_tank_config.state = new_pressurant_tank_state

        # update the time
        time_s += dt_s

        # record the state at this timestep
        sim_points.append(SimPoint(
            time_s=time_s,
            engine=None,
            tanks={
                'liquid_tank': liquid_tank_config.state,
                'pressurant_tank': pressurant_tank_config.state
            },
            injectors_mdot={'injector': liquid_mdot_kg_s},
            regulators_mdot={'regulator': pressurant_mdot_kg_s}
        ))


        # end simulation if liquid dryout occurs
        if liquid_dryout:
            print(f"Liquid dryout occurred at time {time_s:.2f} s, ending simulation.")
            break

        # print data for debugging
        if liquid_tank_config.state.pressure_pa is not None:
            print(f"time: {time_s:.2f} s,   liquid tank pressure: {liquid_tank_config.state.pressure_pa/1e5:.2f} bar,  liquid mass: {liquid_tank_config.state.liquid_mass_kg:.2f} kg,    liquid mass out: {liquid_mdot_kg_s*dt_s:.3f} kg,    pressurant mass out: {pressurant_mdot_kg_s*dt_s:.3f} kg")
        
    if record:
        return SimRecord(points=sim_points)
    else:
        return None
    
