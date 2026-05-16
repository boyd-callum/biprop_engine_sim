# Biprop Engine Sim

A Python simulation for a pressure-fed bipropellant rocket engine. The sim models tank blowdown, injector flow, regulated propellant pressurisation, chamber combustion properties using RocketCEA, and thrust over time.

The current focus is higher-level engine and system behaviour instead of a high-fidelity transient combustion simulator. It is mainly intended for early sizing, designing feed systems, and understanding how the state of the tanks and engine change during the burn.

The code is set up around generic fluid objects. It has mainly been tested with nitrous oxide, ethanol, and nitrogen so far, but other fluids can be added by defining a new `Fluid` object and selecting it in a simulation case.

## Current Features

- General fluid/tank framework using CoolProp fluid objects
- Fluid definitions can be swapped by adding a new `Fluid` object and specifying it in a case
- Tested so far with:
  - nitrous oxide as a self-pressurised oxidiser
  - ethanol as a pressurised liquid fuel
  - nitrogen as a pressurant gas
- Self-pressurised tank model for saturated/blowdown fluids
- Pressurised liquid tank model with separate pressurant gas tracking
- Gas bottle model for pressurant storage
- Regulator mass-flow demand model
- Liquid, gas, and two-phase injector flow models
- RocketCEA-based chamber and nozzle performance
- Dynamic plotting of tank, injector, regulator, and engine outputs
- CSV logging for post-processing

## Model Overview

At each timestep, the simulation:

1. Advances the oxidiser tank and injector flow
2. Advances the fuel tank and injector flow
3. Calculates pressurant demand from any pressurised liquid tanks
4. Advances the pressurant tank
5. Calculates chamber state from propellant mass flow rates
6. Records tank, injector, regulator, and engine outputs

The chamber model currently uses a quasi-steady relationship:

```text
Pc = mdot * cstar / At
```

where `cstar` is calculated using RocketCEA and corrected using a configurable c-star efficiency.

Thrust is then calculated using:

```text
F = Cf * Pc * At
```

where `Cf` is also calculated using RocketCEA and corrected using a configurable thrust coefficient efficiency.

## Project Structure

```text
biprop_engine_sim/
├── main.py               # Entry point for running and plotting a case
├── cases.py              # Simulation case definitions
├── simulation.py         # Main biprop simulation loop
├── timestep_advance.py   # Tank timestep update functions
├── tanks.py              # Tank configs, states, and state solvers
├── injector.py           # Injector config and mass-flow models
├── regulator.py          # Regulator config and pressurant demand model
├── engine.py             # RocketCEA engine/chamber/nozzle model
├── fluid.py              # CoolProp fluid wrapper
├── outputs.py            # Plotting and CSV logging helpers
├── constants.py          # Shared constants
└── helpers.py            # General helper functions
```

## Setup

Create and activate a virtual environment. Tested only using Python 3.14.3.

```powershell
python -m venv .venv
.venv\Scripts\activate
```

Install dependencies:

```powershell
python -m pip install -r requirements.txt
```


## Running a Simulation

Run:

```powershell
python main.py
```

The default case is selected in `main.py`, usually from `cases.py`:

```python
case = full_biprop_case
sim_record = biprop_simulate(case, record=True)
```

After the run, the script can:

- print live timestep status
- save plots
- save a CSV file of simulation results

Example output:

```text
t=0.010 s | Pc=42.50 bar | F=3142.0 N | N2O=59.95 bar | EtOH=50.00 bar | N2=299.44 bar
```

## Defining a Case

Simulation cases are defined in `cases.py`.

A full biprop case generally needs:

- one oxidiser tank
- one fuel tank
- optionally one pressurant tank
- one oxidiser injector
- one fuel injector
- optionally one regulator
- one engine config
- simulation settings




Depending on the propellants in the tanks, they can be configured differently. 

```python
n2o_tank = TankConfig(
    name="Nitrous Oxide Tank",
    role="oxidiser",
    fluid=NITROUS_OXIDE,
    tank_volume_m3=0.0088,    # 8.8L
    phase_model="self_pressurised",
)
```

Each tank also needs an initial conditon. These can be set through various modes:

- Pressure and Mass
  - Good for self-pressurised propellants, like N2O
- Pressure and Temperature
  - Good for single-phase (generally gas) propellants
- Pressure, Temperature, and Mass
  - Needed for liquid tanks with a seperate pressurant gas


```python
n2o_tank_initial = TankInitialCondition(
    mode="pressure_mass",
    pressure_pa=60e+5,    # 60 bar
    total_mass_kg=5.85,
)  
```

The fuel and oxidiser injectors also need their roles set correctly:

```python
n2o_injector = InjectorConfig(
    cd=0.78,
    area_m2=1.05e-5,
    role="oxidiser",
)

ethanol_injector = InjectorConfig(
    cd=0.72,
    area_m2=3.20e-6,
    role="fuel",
)
```

The engine config links the feed system to RocketCEA:

```python
engine_config = EngineConfig(
    geometry=EngineGeometry(
        nozzle_throat_area_m2=2.20e-4,
        expansion_ratio=4.5,
    ),
    cstar_efficiency=0.94,
    cf_efficiency=0.97,
)
```

## Outputs

The plotting and logging helpers generate time histories for:

- tank pressure
- tank temperature
- liquid mass
- vapour mass
- pressurant gas mass
- total tank mass
- internal energy
- injector mass flow rates
- regulator mass flow rates
- chamber pressure
- chamber temperature
- thrust
- specific impulse
- mixture ratio
- c-star
- thrust coefficient

CSV results can be loaded into MATLAB, Excel, Python, or other analysis tools.

## Current Limitations

This is still a development model. Some known limitations:

- Chamber pressure is currently calculated explicitly from injector mass flow, rather than by solving the coupled injector/nozzle pressure balance implicitly.
- Nitrous oxide dryout can cause numerical issues if the model switches sharply between liquid/two-phase and gas-only injector flow.
- The regulator model is simplified and does not yet fully enforce realistic upstream/downstream pressure limits.
- The pressurised liquid tank model is sensitive to state solver tolerances.
- Combustion is quasi-steady - chamber filling, ignition transients, combustion instability, and real combustion speeds are not modelled.
- Heat transfer to the chamber wall, tank walls, and feed lines is not currently included.


## ToDo:

- Improve regulator model to include flow curves from manufacturers
- Implement sweeps and goal seeking for better engine sizing
- Validate against actual engines with various different propellants
