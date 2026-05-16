
from dataclasses import dataclass
from typing import Literal


from injector import InjectorConfig
from tanks import TankState, TankConfig
from engine import EngineState, EngineConfig
from cases import SimCase
from regulator import RegulatorConfig
from constants import *
from helpers import get_single_config_by_role, format_bar, format_force, get_tank_debug_string
from timestep_advance import advance_propellant_tank, advance_pressurant_tank


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





# --------------------
# Full simulations

# def single_tank_simulate(
#     case: SimCase,
#     record: bool = True
# ) -> SimRecord | None:
    
#     """
#     Generic simulator for a single tank blowdown case, with a single injector
#     """

#     if case.tank_configs is None or len(case.tank_configs) == 0:
#         raise ValueError(f"No tank configs included in case {case.name}.")
#     if case.injector_configs is None or len(case.injector_configs) == 0:
#         raise ValueError(f"No injector configs included in case {case.name}.")
    
#     if len(case.tank_configs) > 1:
#         raise ValueError(f"single_tank_simulate is only for simulating one tank. {len(case.tank_configs)} tank configs where given.")
#     if len(case.injector_configs) > 1:
#         raise ValueError(f"single_tank_simulate is only for simulating one injector. {len(case.injector_configs)} injector configs where given.")
#     if len(case.tank_initial_conditions) > 1:
#         raise ValueError(f"single_tank_simulate is only for simulating one tank. {len(case.injector_configs)} tank initial conditions where given.")


#     # get the first entry in both tank and injector congifs. as single tank / injector, should be only entry
#     tank_config = list(case.tank_configs.values())[0]
#     injector_config = list(case.injector_configs.values())[0]
#     initial_condition = list(case.tank_initial_conditions.values())[0]

#     # initialise the tank state from the initial condition
#     tank_config.state = tank_config.initialise_tank_state(initial_condition)

#     dt_s = case.settings.dt_s
#     t_final_s = case.settings.t_final_s

#     time_s = 0.0
#     sim_points = []

#     # save the initial state to the record
#     sim_points.append(SimPoint(
#         time_s = time_s,
#         engine = None,
#         tanks = {'tank' : tank_config.state},
#         injectors_mdot={'injector' : 0.0}
#     ))

#     # While loop to advance the simulation in time
#     while time_s < t_final_s:
        
#         # Advance the simulation by one timestep
#         try:

#             if tank_config.phase_model == "self_pressurised":
#                 new_tank_state = blowdown_advance_timestep(
#                     tank_config=tank_config,
#                     injector_config=injector_config,
#                     dt_s=dt_s
#                 )
#             elif tank_config.phase_model == "single_phase":
#                 new_tank_state = gas_advance_timestep(
#                     tank_config=tank_config,
#                     injector_config=injector_config,
#                     dt_s=dt_s
#                 )
#             else:
#                 raise NotImplementedError(f"Phase model {tank_config.phase_model} has not yet been implmented.")

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
#             tanks={'tank': tank_config.state},
#             injectors_mdot={'injector': injector_mdot_kg_s}
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
    




# def pressurised_liquid_simulate(
#         case: SimCase,
#         record: bool = True
# ) -> SimRecord | None:
#     """
#     Simulator for a pressurised liquid tank with a single injector, using the given case.
#     """
    
#     if case.tank_configs is None or len(case.tank_configs) != 2:
#         raise ValueError(f"Pressurised liquid simulate requires exactly 2 tank configs, one for the liquid tank and one for the pressurant tank. {len(case.tank_configs) if case.tank_configs is not None else 0} were given.")
#     if case.injector_configs is None or len(case.injector_configs) != 1:
#         raise ValueError(f"Pressurised liquid simulate requires exactly 1 injector config. {len(case.injector_configs) if case.injector_configs is not None else 0} were given.")
#     if case.regulator_configs is None or len(case.regulator_configs) != 1:
#         raise ValueError(f"Pressurised liquid simulate requires exactly 1 regulator config. {len(case.regulator_configs) if case.regulator_configs is not None else 0} were given.")
#     if case.tank_initial_conditions is None or len(case.tank_initial_conditions) != 2:
#         raise ValueError(f"Pressurised liquid simulate requires exactly 2 tank initial conditions, one for the liquid tank and one for the pressurant tank. {len(case.tank_initial_conditions) if case.tank_initial_conditions is not None else 0} were given.")
    

#     # determine which tank config is the liquid and which is the pressurant based on the "role" field in the tank config
#     liquid_tank_config = None
#     pressurant_tank_config = None
#     for tank_name, tank_config in case.tank_configs.items():
#         if tank_config.role == "fuel" or tank_config.role == "oxidiser":
#             liquid_tank_config = tank_config
#         elif tank_config.role == "pressurant":
#             pressurant_tank_config = tank_config
#         else:
#             raise ValueError(f"Tank config {tank_name} has invalid role {tank_config.role}. Role must be either 'liquid' or 'pressurant'.")
        
    
#     if liquid_tank_config is None:
#         raise ValueError("No liquid tank config found in case.")
#     if pressurant_tank_config is None:
#         raise ValueError("No pressurant tank config found in case.")
    
#     # determine which initial condition corresponds to which tank
#     liquid_tank_initial_condition = None
#     pressurant_tank_initial_condition = None
#     for tank_name, initial_condition in case.tank_initial_conditions.items():
#         # print(f"tank name: {tank_name}, liquid tank name: {liquid_tank_config.name}, pressurant tank name: {pressurant_tank_config.name}")
#         if tank_name == liquid_tank_config.name:
#             liquid_tank_initial_condition = initial_condition
#         elif tank_name == pressurant_tank_config.name:
#             pressurant_tank_initial_condition = initial_condition
#         else:
#             raise ValueError(f"Tank initial condition {tank_name} does not match any tank config names.")
    
#     if liquid_tank_initial_condition is None:
#         raise ValueError("No initial condition found for liquid tank.")
#     if pressurant_tank_initial_condition is None:
#         raise ValueError("No initial condition found for pressurant tank.")


#     injector_config = list(case.injector_configs.values())[0]
#     regulator_config = list(case.regulator_configs.values())[0]


#     # initialise the tank states from the initial conditions
#     pressurant_tank_config.state = pressurant_tank_config.initialise_tank_state(pressurant_tank_initial_condition)
#     liquid_tank_config.state = liquid_tank_config.initialise_tank_state(liquid_tank_initial_condition)

#     # sim settings
#     dt_s = case.settings.dt_s
#     t_final_s = case.settings.t_final_s

#     time_s = 0.0
#     sim_points = []

#     # record the initial state
#     sim_points.append(SimPoint(
#         time_s = time_s,
#         engine = None,
#         tanks = {
#             'liquid_tank' : liquid_tank_config.state,
#             'pressurant_tank' : pressurant_tank_config.state
#         },
#         injectors_mdot={'injector' : 0.0}
#     ))


#     # main while loop to advance the simulation in time
#     while time_s < t_final_s:

#         # advance the simulation by one timestep
#         try:

#             (
#                 new_liquid_tank_state,
#                 new_pressurant_tank_state,
#                 liquid_mdot_kg_s,
#                 pressurant_mdot_kg_s,
#                 liquid_dryout
#             ) = pressurised_liquid_advance_timestep(
#                 liquid_tank_config=liquid_tank_config,
#                 injector_config=injector_config,
#                 regulator_config=regulator_config,
#                 pressurant_tank_config=pressurant_tank_config,
#                 dt_s=dt_s
#             )
        
#         except RuntimeError as e:
#             print(f"Error advancing timestep at time {time_s:.2f} s: {e}")
#             print(f"final liquid tank state at time of error: pressure={liquid_tank_config.state.pressure_pa:.2f} Pa, temperature={liquid_tank_config.state.temperature_k:.2f} K, total_mass={liquid_tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={liquid_tank_config.state.total_internal_energy_j:.2f} J, liquid_mass={liquid_tank_config.state.liquid_mass_kg:.2f} kg, pressurant_mass={liquid_tank_config.state.pressurant_gas_mass_kg:.2f} kg")
#             print(f"final pressurant tank state at time of error: pressure={pressurant_tank_config.state.pressure_pa:.2f} Pa, temperature={pressurant_tank_config.state.temperature_k:.2f} K, total_mass={pressurant_tank_config.state.total_mass_kg:.2f} kg, total_internal_energy={pressurant_tank_config.state.total_internal_energy_j:.2f} J")
#             break


        
#         # update the tank states
#         liquid_tank_config.state = new_liquid_tank_state
#         pressurant_tank_config.state = new_pressurant_tank_state

#         # update the time
#         time_s += dt_s

#         # record the state at this timestep
#         sim_points.append(SimPoint(
#             time_s=time_s,
#             engine=None,
#             tanks={
#                 'liquid_tank': liquid_tank_config.state,
#                 'pressurant_tank': pressurant_tank_config.state
#             },
#             injectors_mdot={'injector': liquid_mdot_kg_s},
#             regulators_mdot={'regulator': pressurant_mdot_kg_s}
#         ))


#         # end simulation if liquid dryout occurs
#         if liquid_dryout:
#             print(f"Liquid dryout occurred at time {time_s:.2f} s, ending simulation.")
#             break

#         # print data for debugging
#         if liquid_tank_config.state.pressure_pa is not None:
#             print(f"time: {time_s:.2f} s,   liquid tank pressure: {liquid_tank_config.state.pressure_pa/1e5:.2f} bar,  liquid mass: {liquid_tank_config.state.liquid_mass_kg:.2f} kg,    liquid mass out: {liquid_mdot_kg_s*dt_s:.3f} kg,    pressurant mass out: {pressurant_mdot_kg_s*dt_s:.3f} kg")
        
#     if record:
#         return SimRecord(points=sim_points)
#     else:
#         return None
    








def biprop_simulate(   
        case: SimCase,
        record: bool = True 
    ) -> SimRecord | None:


    if case.tank_configs is None:
        raise ValueError("biprop_simulate requires tanks.")
    if case.injector_configs is None:
        raise ValueError("biprop_simulate requires injectors.")
    if case.tank_initial_conditions is None:
        raise ValueError("biprop_simulate requires tank initial conditions.")
    if case.engine_config is None:
        raise ValueError("biprop_simulate requires an engine config.")
    if case.settings is None:
        raise ValueError("biprop_simulate requires simulation settings.")
    

    # -----------------------
    # 1. get required tanks by role

    fuel_tank_id, fuel_tank_config = get_single_config_by_role(
        configs=case.tank_configs,
        role="fuel",
        config_type_name="tank"
    )
    

    ox_tank_id, ox_tank_config = get_single_config_by_role(
        configs=case.tank_configs,
        role="oxidiser",
        config_type_name="tank"
    )

    # ------------------
    # 2. get required injectors by role

    fuel_injector_id, fuel_injector_config = get_single_config_by_role(
        configs=case.injector_configs,
        role="fuel",
        config_type_name="injector"
    )

    ox_injector_id, ox_injector_config = get_single_config_by_role(
        configs=case.injector_configs,
        role="oxidiser",
        config_type_name="injector"
    )

    # --------------------
    # 3. determine if we need pressurant tank

    propellant_tank_configs = {
        fuel_tank_id: fuel_tank_config,
        ox_tank_id: ox_tank_config
    }

    pressurised_propellant_tank_ids = [
        tank_id 
        for tank_id, tank_confifg in propellant_tank_configs.items()
        if tank_confifg.phase_model == "pressurised_liquid"
    ]

    pressurant_tank_id = None
    pressurant_tank_config = None
    regulator_id = None
    regulator_config = None

    if len(pressurised_propellant_tank_ids) > 0:

        pressurant_tank_id, pressurant_tank_config = get_single_config_by_role(
            configs=case.tank_configs,
            role="pressurant",
            config_type_name="tank",
        )

        if case.regulator_configs is None:
            raise ValueError(
                "At least one propellant tank is pressurised_liquid, "
                "but no regulator configs were provided."
            )

        if len(pressurised_propellant_tank_ids) > 1:
            raise NotImplementedError(
                "Multiple pressurised propellant tanks are present. "
                "This needs an explicit regulator-to-tank mapping."
            )

        if len(case.regulator_configs) != 1:
            raise ValueError(
                "This version expects exactly one regulator config when one "
                "propellant tank is pressurised_liquid."
            )

        regulator_id, regulator_config = next(iter(case.regulator_configs.items()))        



    # -------------------------
    # 4. initialise tank states


    tank_states: dict[str, TankState] = {}

    for tank_id, tank_config in case.tank_configs.items():

        initial_condition = case.tank_initial_conditions.get(tank_id)
        
        if initial_condition is None:
            raise ValueError(
                f"No initial condition found for tank '{tank_id}' "
                f"with name '{tank_config.name}'."
            )
    

        tank_state = tank_config.initialise_tank_state(initial_condition)

        tank_config.state = tank_state
        tank_states[tank_id] = tank_state

    

    #------------------
    # 5. initialise engine state

    engine = case.engine_config

    engine.fuel_cea_name = fuel_tank_config.fluid.cea_name
    engine.ox_cea_name = ox_tank_config.fluid.cea_name

    if engine.state is None:
        raise ValueError("Engine is missing state.")

    if engine.state.chamber_pressure_pa is None:
        raise ValueError("Engine state is missing pressure after initialisation.")

    # -------------------
    # 6. set up main loop


    sim_points: list[SimPoint] = []

    time_s = 0.0
    dt_s = case.settings.dt_s
    end_time_s = case.settings.t_final_s

    # for printing live data
    step_index = 0


    # -----------
    # 7. main loop

    while time_s <= end_time_s:

        tank_updates: dict[str, TankState] = {}
        injector_mdots: dict[str, float] = {}
        regulator_mdots: dict[str, float] = {}

        total_pressurant_mdot_kg_s = 0.0

        downstream_pressure_pa = engine.state.chamber_pressure_pa


        # -------------
        # 7a. advance ox path

        ox_uses_pressurisation = ox_tank_config.phase_model == "pressurised_liquid"

        ox_result = advance_propellant_tank(
            tank_config=ox_tank_config,
            injector_config=ox_injector_config,
            dt_s=dt_s,
            downstream_pressure_pa=downstream_pressure_pa,
            pressurant_tank_config=pressurant_tank_config if ox_uses_pressurisation else None,
            regulator_config= regulator_config if ox_uses_pressurisation else None
        )


        tank_updates[ox_tank_id] = ox_result.tank_update

        if ox_result.injector_mdot_kg_s is not None:
            injector_mdots[ox_injector_id] = ox_result.injector_mdot_kg_s
        

        if ox_result.regulator_mdot_kg_s is not None:
            
            if regulator_id is None:
                raise ValueError("Oxidiser regulator mdot exists, but regulator_id is None.")
            
            regulator_mdots[regulator_id] = ox_result.regulator_mdot_kg_s
            total_pressurant_mdot_kg_s += ox_result.regulator_mdot_kg_s


        #---------------------------
        # 7b. advance fuel path


        fuel_uses_pressurisation = fuel_tank_config.phase_model == "pressurised_liquid"

        fuel_result = advance_propellant_tank(
            tank_config=fuel_tank_config,
            injector_config=fuel_injector_config,
            dt_s=dt_s,
            downstream_pressure_pa=downstream_pressure_pa,
            pressurant_tank_config=pressurant_tank_config if fuel_uses_pressurisation else None,
            regulator_config=regulator_config if fuel_uses_pressurisation else None
        )


        tank_updates[fuel_tank_id] = fuel_result.tank_update

        if fuel_result.injector_mdot_kg_s is not None:
            injector_mdots[fuel_injector_id] = fuel_result.injector_mdot_kg_s
        

        if fuel_result.regulator_mdot_kg_s is not None:

            if regulator_id is None:
                raise ValueError("Fuel regulator mdot exists, but regulator_id is None.")
            
            regulator_mdots[regulator_id] = fuel_result.regulator_mdot_kg_s
            total_pressurant_mdot_kg_s += fuel_result.regulator_mdot_kg_s

            
        

        # -------------
        # 7c. advance pressurant tank, only once both prop tanks are advanced

        if total_pressurant_mdot_kg_s > 0.0:

            if pressurant_tank_id is None or pressurant_tank_config is None:
                raise ValueError(
                    "Pressurant mdot is nonzero, but no pressurant tank is available."
                )
            
            new_pressurant_tank_state = advance_pressurant_tank(
                tank_config=pressurant_tank_config,
                dt_s=dt_s,
                pressurant_mdot_kg_s=total_pressurant_mdot_kg_s
            )


            tank_updates[pressurant_tank_id] = new_pressurant_tank_state
        

        #-----------------------
        # 7d. apply tank updates

        for tank_id, new_tank_state in tank_updates.items():

            case.tank_configs[tank_id].state = new_tank_state
            tank_states[tank_id] = new_tank_state
        

        # ------------------
        # 7e. calculate engine state

        try:
            ox_mdot_kg_s = injector_mdots.get(ox_injector_id)
            fuel_mdot_kg_s = injector_mdots.get(fuel_injector_id)

            if ox_mdot_kg_s is None:
                raise ValueError(
                    f"No oxidiser mdot found for injector '{ox_injector_id}'."
                )

            if fuel_mdot_kg_s is None:
                raise ValueError(
                    f"No fuel mdot found for injector '{fuel_injector_id}'."
                )

            engine.state = engine.calculate_state(
                ox_mdot_kg_s=ox_mdot_kg_s,
                fuel_mdot_kg_s=fuel_mdot_kg_s,
                ambient_pressure_pa=ATMOSPHERE_PRESSURE_PA,
                previous_state=engine.state,
            )

        except Exception as exc:
            raise RuntimeError(
                "Failed while calculating engine state. "
                f"step={step_index}, t={time_s:.3f} s, "
                f"ox_mdot={injector_mdots.get(ox_injector_id)}, "
                f"fuel_mdot={injector_mdots.get(fuel_injector_id)}, "
                f"previous_pc={format_bar(downstream_pressure_pa)}"
            ) from exc

        # ----------
        # 7f. record timestep


        time_s += dt_s

        if record:

            sim_points.append(
                SimPoint(
                    time_s=time_s,
                    engine=engine.state,
                    tanks=tank_states.copy(),
                    injectors_mdot=injector_mdots.copy(),
                    regulators_mdot=regulator_mdots.copy()
                )
            )


        ox_tank_state = tank_states[ox_tank_id]
        fuel_tank_state = tank_states[fuel_tank_id]

        if step_index % case.settings.print_steps == 0:

        
            print_parts = [
                f"t={time_s:.3f} s",
                f"Pc={format_bar(engine.state.chamber_pressure_pa)} bar",
                f"F={format_force(engine.state.thrust_n)} N",
                f"N2O={format_bar(ox_tank_state.pressure_pa)} bar",
                f"EtOH={format_bar(fuel_tank_state.pressure_pa)} bar",
            ]

            if pressurant_tank_id is not None:
                pressurant_tank_state = tank_states[pressurant_tank_id]
                print_parts.append(
                    f"N2={format_bar(pressurant_tank_state.pressure_pa)} bar"
                )

            print(" | ".join(print_parts))


        min_pressure = ATMOSPHERE_PRESSURE_PA*1.5

        if ox_tank_state.pressure_pa is None or fuel_tank_state.pressure_pa is None:
            raise ValueError("Ox or Fuel pressure is None")

        if ox_tank_state.pressure_pa < min_pressure:
            print(f"Ox tank pressure dropped below {format_bar(min_pressure)} bar, sim ended at {time_s} seconds.")
            break

        if fuel_tank_state.pressure_pa < min_pressure:
            print(f"Fuel tank pressure dropped below {format_bar(min_pressure)} bar, sim ended at {time_s} seconds.")
            break
            
        if engine.state.thrust_n == 0.0:
            print(f"Engine stopped producing thrust, sim ended at {time_s} seconds.")
            break
        
        step_index += 1
        
        

    if record:
        return SimRecord(
            points = sim_points
        )
    
    return None