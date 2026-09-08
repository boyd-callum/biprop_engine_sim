from __future__ import annotations


from dataclasses import dataclass
from typing import Literal, Callable

from constants import *
from fluid import Fluid
from helpers import bisection_search, brent_search

# -------------------------------
# Solver settings

# the reconstructed vapour must fill the ullage to within this
SATURATED_VOLUME_TOLERANCE_M3 = 1e-7

# initial half-width of the temperature bracket, as a fraction of the guess
TEMPERATURE_BRACKET_FRACTION = 0.05

# how far the bracket is allowed to widen before giving up
TEMPERATURE_BRACKET_MAX_EXPANSIONS = 12


def solve_temperature_at_density(
    fluid: Fluid,
    residual_func: Callable[[float], float],
    temperature_guess_k: float,
    tank_name: str = "tank",
) -> float:
    """
    Solve for the temperature that zeroes a residual, bracketing outward from a
    guess and never evaluating outside the fluid model's valid range.

    CoolProp extrapolates silently outside that range instead of raising, so a
    solver given wide fixed bounds can converge on physically meaningless states.
    Bounds here come from the fluid's own limits and the bracket starts narrow
    around the previous temperature, which is normally within a degree or two.
    """

    minimum_temperature_k = fluid.get_Tmin()
    maximum_temperature_k = fluid.get_Tmax()

    temperature_guess_k = min(
        max(temperature_guess_k, minimum_temperature_k),
        maximum_temperature_k
    )

    half_width_k = max(temperature_guess_k * TEMPERATURE_BRACKET_FRACTION, 1.0)

    low_k = temperature_guess_k
    high_k = temperature_guess_k

    for _ in range(TEMPERATURE_BRACKET_MAX_EXPANSIONS):

        low_k = max(temperature_guess_k - half_width_k, minimum_temperature_k)
        high_k = min(temperature_guess_k + half_width_k, maximum_temperature_k)

        if residual_func(low_k) * residual_func(high_k) <= 0.0:
            return brent_search(
                residual_func=residual_func,
                bounds=[low_k, high_k],
                tolerance=1e-3,
            )

        if low_k <= minimum_temperature_k and high_k >= maximum_temperature_k:
            break

        half_width_k *= 2.0

    raise ValueError(
        f"Could not bracket a temperature for {tank_name} between "
        f"{minimum_temperature_k:.2f} K and {maximum_temperature_k:.2f} K "
        f"(guess was {temperature_guess_k:.2f} K). The state is probably outside "
        f"the range of the fluid model."
    )


# -------------------------------
# Dataclasses for the tanks

SourceRole = Literal["fuel", "oxidiser", "pressurant"]
PhaseModel = Literal["pressurised_liquid", "single_phase", "self_pressurised", "unknown"]
TankInitMode = Literal[
    "pressure_mass",         # good for self-pressurised N2O
    "pressure_temperature",  # good for simple single-phase (gas) tanks
    "pressure_temperature_mass",      # good for pressurised liquid tanks
]

@dataclass
class TankInitialCondition:
    mode: TankInitMode
    pressure_pa: float | None = None
    total_mass_kg: float | None = None
    temperature_k: float | None = None

    def __post_init__(self):
        if self.mode == "pressure_mass":
            if self.pressure_pa is None or self.total_mass_kg is None:
                raise ValueError("pressure_mass requires pressure_pa and total_mass_kg")
        
        elif self.mode == "pressure_temperature":
            if self.pressure_pa is None or self.temperature_k is None:
                raise ValueError("pressure_temperature requires pressure_pa and temperature_k")

        elif self.mode == "pressure_temperature_mass":
            if self.pressure_pa is None or self.temperature_k is None or self.total_mass_kg is None:
                raise ValueError("pressure_temperature_mass requires pressure_pa, temperature_k, and total_mass_kg")

@dataclass
class TankConfig:
    name: str               # eg "ethanol tank", "n2o tank"
    role: SourceRole        # fuel / oxidiser / pressurant
    fluid: Fluid
    tank_volume_m3: float
    phase_model: PhaseModel = 'unknown'
    state: TankState | None = None
    pressurant_fluid: Fluid | None = None


    def get_saturated_mass_split(
            self,
            total_mass_kg: float,
            saturation_properties: dict
    ) -> tuple[float, float]:
        """
        Equilibrium liquid/vapour mass split for a given total mass at a given
        saturation state.

        The raw split goes non-physical outside the dome - negative liquid mass
        above the saturated-vapour line, negative vapour mass below the saturated
        liquid line - so it is clamped to the physical range. That stops the energy
        calculation being built on a negative mass.

        Note the clamped energy is NOT monotonic in temperature: once the split
        pins to all-vapour the energy follows ug(T), which turns over near the
        critical point. A root found in the clamped region is therefore meaningless,
        and callers must validate the result - the volume check in
        _state_from_mass_and_energy_self_pressurised is what rejects it, after
        which reconstruct_blowdown_state falls back to the single-phase solver.
        """

        vf = saturation_properties['vf']
        vg = saturation_properties['vg']

        mass_vapour_kg = (self.tank_volume_m3 - total_mass_kg * vf) / (vg - vf)

        mass_vapour_kg = min(max(mass_vapour_kg, 0.0), total_mass_kg)
        mass_liquid_kg = total_mass_kg - mass_vapour_kg

        return mass_liquid_kg, mass_vapour_kg


    def calculate_internal_energy(
            self,
            total_mass_kg: float,
            temperature_k: float
    ) -> float:
        """
        Find the total internal energy of a two-phase saturated state with the specified total mass and temperature.
        """

        # find the saturation properties at the given temperature
        saturation_properties = self.fluid.get_saturation_properties_from_temp(temperature_k)

        mass_liquid_kg, mass_vapour_kg = self.get_saturated_mass_split(
            total_mass_kg=total_mass_kg,
            saturation_properties=saturation_properties
        )

        # calculate the internal energy at the given temperature
        total_internal_energy_j = (
            mass_liquid_kg * saturation_properties['uf']
            + mass_vapour_kg * saturation_properties['ug']
        )


        return total_internal_energy_j
    
# -----------------------------
# Functions to initialise the state of the tanks before the first timestep

    def _initialise_self_pressurised_tank_from_pressure_mass(
            self, 
            initial_condition: TankInitialCondition
    ) -> TankState:
        
        """
        Initialise a saturated self-pressurised tank from total mass and tank pressure.

        Assumes the tank is in liquid-vapour equilibrium and that the specified
        pressure lies between the fluid triple and critical pressures.
        """


        if initial_condition.pressure_pa is None or initial_condition.total_mass_kg is None:
            raise ValueError("pressure_mass requires pressure_pa and total_mass_kg")

        pressure_pa = initial_condition.pressure_pa
        total_mass_kg = initial_condition.total_mass_kg

        # find the saturation temperature given the inital pressure
        saturation_properties = self.fluid.get_saturation_properties_from_pressure(pressure_pa)
        saturation_temp_k = saturation_properties['tsat']

        vf = saturation_properties['vf']    # liquid specific volume (m3 / kg)
        vg = saturation_properties['vg']    # vapour specific volume (m3 / kg)


        mass_vapour_kg = (self.tank_volume_m3 - total_mass_kg*vf)/ (vg - vf)
        mass_liquid_kg = total_mass_kg - mass_vapour_kg

        # validate the phase split
        if mass_vapour_kg < 0.0 or mass_vapour_kg > total_mass_kg:
            raise ValueError(f"Initial condition for {self.name} is not consistent with a two-phase saturated state at {pressure_pa} Pa.")

        liquid_volume_m3 = mass_liquid_kg * vf
        ullage_volume_m3 = self.tank_volume_m3 - liquid_volume_m3   
        
        # validate the volumes
        vapour_volume_m3 = mass_vapour_kg * vg
        volume_error = abs(vapour_volume_m3 - ullage_volume_m3)
        if volume_error > 1e-7:
            raise ValueError(f"Volume calculation failed. Volume error of {volume_error} m^3.")

        # calculate the total internal energy
        liquid_energy_j = mass_liquid_kg * saturation_properties['uf']
        vapour_energy_j = mass_vapour_kg * saturation_properties['ug']
        total_internal_energy_j = liquid_energy_j + vapour_energy_j



        return TankState(
            config=self,
            pressure_pa=pressure_pa,
            temperature_k=saturation_temp_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=mass_liquid_kg,
            vapour_mass_kg=mass_vapour_kg,
            pressurant_gas_mass_kg=None,
            ullage_volume_m3=ullage_volume_m3,
            liquid_volume_m3=liquid_volume_m3
        )


    def _initialise_single_phase_tank_from_pressure_temperature(
            self,
            initial_condition: TankInitialCondition
    ) -> TankState:
        """
        Initialise a gas tank from specified pressure and temperature, using real gas behaviour
        """

        
        if initial_condition.pressure_pa is None or initial_condition.temperature_k is None:
            raise ValueError("pressure_temperature requires pressure_pa and temperature_k")

        pressure_pa = initial_condition.pressure_pa
        temperature_k = initial_condition.temperature_k

        # Get the real fluid properties at current state
        density_kg_m3 = self.fluid.get_fluid_density_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )
        specific_internal_energy_j_kg = self.fluid.get_specific_internal_energy_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )

        # Convert intensive properties into extensive tank quantities.
        total_mass_kg = density_kg_m3 * self.tank_volume_m3
        total_internal_energy_j = specific_internal_energy_j_kg * total_mass_kg

        

        return TankState(
            config=self,
            pressure_pa=initial_condition.pressure_pa,
            temperature_k=initial_condition.temperature_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=0.0,
            vapour_mass_kg=0.0,
            pressurant_gas_mass_kg=total_mass_kg,
            ullage_volume_m3=self.tank_volume_m3,
            liquid_volume_m3=0.0
        )


    def _initialise_liquid_tank_from_pressure_temperature_mass(
            self,
            initial_condition: TankInitialCondition
    ) -> TankState:
        
        """
        initialise a gas-pressurised liquid tank from:
            - tnak pressure
            - tank temperature
            - liquid mass

        assumptions:
            - liquid and ullage gas are initially in thermal equilibirium
            - ullage gas fills all volume not occupied by liquid
            - pressurant gas is single-phase (using real-fluid properties)
        """


        if initial_condition.pressure_pa is None:
            raise ValueError("Initial condition is missing pressure_pa.")
        if initial_condition.temperature_k is None:
            raise ValueError("Initial condition is missing temperature_k.")
        if initial_condition.total_mass_kg is None:
            raise ValueError("Initial condition is missing mass_kg for liquid mass.")

        pressure_pa = initial_condition.pressure_pa
        temperature_k = initial_condition.temperature_k
        liquid_mass_kg = initial_condition.total_mass_kg

        # find liquid density at inital state
        liquid_density_kg_m3 = self.fluid.get_fluid_density_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )

        liquid_volume_m3 = liquid_mass_kg / liquid_density_kg_m3
        ullage_volume_m3 = self.tank_volume_m3 - liquid_volume_m3

        # check to make sure no error in volumes

        if ullage_volume_m3 < 0.0:
            raise ValueError(
                f"Liquid volume exceeds tank volume in {self.name}. "
                f"Liquid volume = {liquid_volume_m3:.6e} m^3, "
                f"tank volume = {self.tank_volume_m3:.6e} m^3"
            )
        
        if self.pressurant_fluid is None:
            raise ValueError(f"Gas-pressurised liquid tank ({self.name}) requires a pressurant fluid.")


        # pressurant gas density and mass in the ullage
        pressurant_density_kg_m3 = self.pressurant_fluid.get_fluid_density_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )

        pressurant_mass_kg = pressurant_density_kg_m3 * ullage_volume_m3

        # specific internal energies
        liquid_u_j_kg = self.fluid.get_specific_internal_energy_from_pressure_temperature(
            T= temperature_k,
            P=pressure_pa
        )
        pressurant_u_j_kg = self.pressurant_fluid.get_specific_internal_energy_from_pressure_temperature(
            T=temperature_k,
            P=pressure_pa
        )


        # find total mass and internal energy
        total_mass_kg = liquid_mass_kg + pressurant_mass_kg
        total_internal_energy_j = liquid_mass_kg*liquid_u_j_kg + pressurant_mass_kg*pressurant_u_j_kg


        return TankState(
            config=self,
            pressure_pa=pressure_pa,
            temperature_k=temperature_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=liquid_mass_kg,
            vapour_mass_kg=0,
            pressurant_gas_mass_kg=pressurant_mass_kg,
            ullage_volume_m3=ullage_volume_m3,
            liquid_volume_m3=liquid_volume_m3,
            pressurant_gas_internal_energy_j=pressurant_mass_kg*pressurant_u_j_kg
        )




    def initialise_tank_state(
        self,
        initial_condition: TankInitialCondition
    ) -> TankState:

        if self.phase_model == "self_pressurised" and initial_condition.mode == "pressure_mass":
            return self._initialise_self_pressurised_tank_from_pressure_mass(initial_condition)

        elif self.phase_model == "single_phase" and initial_condition.mode == "pressure_temperature":
            return self._initialise_single_phase_tank_from_pressure_temperature(initial_condition)

        elif self.phase_model == "pressurised_liquid" and initial_condition.mode == "pressure_temperature_mass":
            return self._initialise_liquid_tank_from_pressure_temperature_mass(initial_condition)

        raise NotImplementedError(f"Tank initialisation not implemented for phase_model={self.phase_model!r}, mode={initial_condition.mode!r}")




#------------------------------------
# Functions to determine thet state of the tank during each timestep after fluid has left


    def _state_from_mass_and_energy_self_pressurised(
            self,
            total_mass_kg: float,
            total_internal_energy_j: float,
            previous_state: TankState | None = None
    ) -> TankState:

        if previous_state is not None and previous_state.temperature_k is not None:
            guessed_temperature_k = previous_state.temperature_k
        else:
            guessed_temperature_k = self.fluid.get_Ttriple() + 10.0   # arbitrary guess above the triple point

        # find the total internal energy at the guessed energy

        guessed_energy_j = self.calculate_internal_energy(total_mass_kg, guessed_temperature_k)

        # now need to root search to find the temperature which gives the correct internal energy for the given mass, using a bisection method

        T_low_k = self.fluid.get_Ttriple()
        T_high_k = self.fluid.get_Tcrit() - 0.01

        energy_error_j = guessed_energy_j - total_internal_energy_j
        iteration = 0

        while abs(energy_error_j) > 1e-3:
            if energy_error_j > 0.0:
                T_high_k = guessed_temperature_k
            else:
                T_low_k = guessed_temperature_k
            
            guessed_temperature_k = 0.5 * (T_low_k + T_high_k)
            guessed_energy_j = self.calculate_internal_energy(total_mass_kg, guessed_temperature_k)
            energy_error_j = guessed_energy_j - total_internal_energy_j
            iteration += 1
        
            if iteration > 100:
                raise RuntimeError(f"Root search failed to converge after 100 iterations. Final energy error of {energy_error_j} J.")
        
        # at this point guessed_temperature_k should be the temperature which gives the correct internal energy for the given mass, so we can calculate the rest of the state properties at this temperature

        correct_temperature_k = guessed_temperature_k

        saturation_properties = self.fluid.get_saturation_properties_from_temp(correct_temperature_k)

        vf = saturation_properties['vf']

        # same clamped split the energy root search used, so the reconstructed
        # state is consistent with the energy that was solved for
        mass_liquid_kg, mass_vapour_kg = self.get_saturated_mass_split(
            total_mass_kg=total_mass_kg,
            saturation_properties=saturation_properties
        )

        liquid_volume_m3 = mass_liquid_kg * vf
        ullage_volume_m3 = self.tank_volume_m3 - liquid_volume_m3

        # The vapour must fill the ullage exactly. If it does not, the solved
        # temperature does not describe a two-phase state at this mass and volume,
        # which means the tank has left the dome - the caller is expected to fall
        # back to the single-phase solver.
        vapour_volume_m3 = mass_vapour_kg * saturation_properties['vg']
        volume_error = abs(vapour_volume_m3 - ullage_volume_m3)

        if volume_error > SATURATED_VOLUME_TOLERANCE_M3:
            raise ValueError(
                f"{self.name} is not a valid two-phase state at "
                f"{correct_temperature_k:.3f} K: vapour volume "
                f"{vapour_volume_m3:.6e} m^3 does not fill the ullage "
                f"{ullage_volume_m3:.6e} m^3 (error {volume_error:.3e} m^3)."
            )


        liquid_energy_j = mass_liquid_kg * saturation_properties['uf']
        vapour_energy_j = mass_vapour_kg * saturation_properties['ug']
        total_internal_energy_j = liquid_energy_j + vapour_energy_j



        return TankState(
            config=self,
            pressure_pa=saturation_properties['psat'],
            temperature_k=correct_temperature_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=mass_liquid_kg,
            vapour_mass_kg=mass_vapour_kg,
            pressurant_gas_mass_kg=None,
            ullage_volume_m3=ullage_volume_m3,
            liquid_volume_m3=liquid_volume_m3
            )
        


    def _state_from_mass_and_energy_single_phase(
            self,
            total_mass_kg: float,
            total_internal_energy_j: float,
            previous_state: TankState | None = None
    ) -> TankState:
        
        """
        find the state of a tank containing single-phase real fluid from its total mass and internal energy. Can be used for:
            - gas
            - liquid
            - supercritical fluid

        Method:
            1. find bulk fluid density from known volume

            2. compute target specific internal energy

            3. solve for temperature that gives this specific internal energy at the computed density

            4. get pressure from this density - temperature state
        
        """

        # check if there is still mass left
        if total_mass_kg <= 0.0:
            raise ValueError(f"mass in {self.name} must remain positive")
        

        # find density
        density_kg_m3 = total_mass_kg / self.tank_volume_m3

        # convert total internal energy into specific internal energy
        target_specific_internal_energy_j_kg = total_internal_energy_j / total_mass_kg

        # use previous temperature as a starting point, as the current temp is probably pretty close
        if previous_state is not None and previous_state.temperature_k is not None:
            temperature_guess_k = previous_state.temperature_k
        else:
            temperature_guess_k = 300 # around ambient

        def temperature_residual(temperature_k: float) -> float:
            """
            residual for the temperature root solve

            at known density, find difference between 
                - internal energy predicted by fluid model at temp
                - target internal energy

            correct temp is the root where this is zero
            """

            specific_internal_energy_j_kg = self.fluid.props_si("U", "D", density_kg_m3, "T", temperature_k)

            residual = specific_internal_energy_j_kg - target_specific_internal_energy_j_kg

            return residual
        

        # Bracket the root around the previous temperature and widen until it
        # straddles, rather than starting from a fixed [10, 1000] K span.
        #
        # The old fixed span was not just wasteful, it was evaluating states the
        # fluid model cannot represent. CoolProp does not raise below its minimum
        # temperature - it silently extrapolates, and returns nonsense
        # (U for nitrogen at D = 300, T = 10 K comes back as -3.1e20 J/kg). The
        # bisection then "worked" only because that garbage happened to sit on the
        # far side of the root.
        temperature_k = solve_temperature_at_density(
            fluid=self.fluid,
            residual_func=temperature_residual,
            temperature_guess_k=temperature_guess_k,
            tank_name=self.name,
        )

        # get pressure from solved density-temp state
        pressure_pa = self.fluid.props_si("P", "D", density_kg_m3, "T", temperature_k)

        
        


        return TankState(
            config=self,
            pressure_pa=pressure_pa,
            temperature_k=temperature_k,
            total_mass_kg=total_mass_kg,
            total_internal_energy_j=total_internal_energy_j,
            liquid_mass_kg=0.0,
            vapour_mass_kg=total_mass_kg,
            pressurant_gas_mass_kg=total_mass_kg if self.role == "pressurant" else None,
            ullage_volume_m3=self.tank_volume_m3,
            liquid_volume_m3=0.0
        )




    def _state_from_mass_and_energy_pressurised_liquid(
            self,
            liquid_mass_kg: float,
            pressurant_gas_mass_kg: float,
            total_internal_energy_j: float,
            previous_state: TankState | None = None
    ) -> TankState:
        


        if self.pressurant_fluid is None:
            raise ValueError(f"{self.name} requires a pressurant fluid.")
        if liquid_mass_kg < 0.0:
            raise ValueError("Liquid mass must be non-negative")
        if pressurant_gas_mass_kg < 0.0:
            raise ValueError("Pressurant gas mass must be non-negative")
        

        if previous_state is not None and previous_state.temperature_k is not None:
            temperature_guess_k = previous_state.temperature_k
        else:
            temperature_guess_k = 300.0



        def solve_pressure_for_temperature(
                temperature_k: float
            ) -> float:

            """
            with a given temperature, finds the pressure for which that temperature is a valid tank state
            """

            def pressure_residual(
                    pressure_pa:float
                )-> float:

                if self.pressurant_fluid is None:
                    raise ValueError

                try:
                    liquid_density_kg_m3 = self.fluid.get_fluid_density_from_pressure_temperature(
                        P=pressure_pa,
                        T=temperature_k
                    )
                    gas_density_kg_m3 = self.pressurant_fluid.get_fluid_density_from_pressure_temperature(
                        P=pressure_pa,
                        T=temperature_k
                    )
                except Exception:
                    return 1e30


                liquid_volume_m3 = liquid_mass_kg / liquid_density_kg_m3
                gas_volume_m3 = pressurant_gas_mass_kg / gas_density_kg_m3

                return self.tank_volume_m3 - liquid_volume_m3 - gas_volume_m3


            pressure_pa = bisection_search(
                residual_func=pressure_residual,
                bounds = [1e4, 1e8],
                tolerance=1e-9
            )
        
            return pressure_pa
        



        
        def temperature_residual(
                temperature_k: float
            )-> float:
            
            if self.pressurant_fluid is None:
                raise ValueError

            pressure_pa = solve_pressure_for_temperature(temperature_k)

            liquid_u_j_kg = self.fluid.get_specific_internal_energy_from_pressure_temperature(
                P=pressure_pa,
                T=temperature_k
            )
            
            gas_u_j_kg = self.pressurant_fluid.get_specific_internal_energy_from_pressure_temperature(
                P=pressure_pa,
                T=temperature_k
            )


            reconstructed_energy_j = liquid_mass_kg*liquid_u_j_kg + pressurant_gas_mass_kg*gas_u_j_kg

            return reconstructed_energy_j - total_internal_energy_j
        
        temp_low_k = max(temperature_guess_k - 40.0, 250.0)
        temp_high_k = min(temperature_guess_k + 40.0, 450.0)

        temperature_k = bisection_search(
            residual_func=temperature_residual,
            bounds=[temp_low_k, temp_high_k],
            tolerance=1e-2
        )

        pressure_pa = solve_pressure_for_temperature(temperature_k)


        liquid_density_kg_m3 = self.fluid.get_fluid_density_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )
        gas_density_kg_m3 = self.pressurant_fluid.get_fluid_density_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )

        liquid_volume_m3 = liquid_mass_kg / liquid_density_kg_m3
        ullage_volume_m3 = self.tank_volume_m3 - liquid_volume_m3

        liquid_u_kg_j = self.fluid.get_specific_internal_energy_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )
        gas_u_kg_j = self.pressurant_fluid.get_specific_internal_energy_from_pressure_temperature(
            P=pressure_pa,
            T=temperature_k
        )

        pressurant_gas_internal_energy_j = pressurant_gas_mass_kg*gas_u_kg_j
        total_internal_energy_j = liquid_mass_kg*liquid_u_kg_j + pressurant_gas_internal_energy_j

        total_mass_kg = liquid_mass_kg + pressurant_gas_mass_kg

        return TankState(
            config = self,
            pressure_pa = pressure_pa,
            temperature_k = temperature_k,
            total_mass_kg = total_mass_kg,
            total_internal_energy_j = total_internal_energy_j,
            liquid_mass_kg = liquid_mass_kg,
            vapour_mass_kg = 0.0,
            pressurant_gas_mass_kg = pressurant_gas_mass_kg,
            ullage_volume_m3 = ullage_volume_m3,
            liquid_volume_m3 = liquid_volume_m3,
            pressurant_gas_internal_energy_j = pressurant_gas_internal_energy_j
        )




    def state_from_mass_and_energy(
        self,
        total_internal_energy_j: float,
        total_mass_kg: float | None = None,
        liquid_mass_kg: float | None = None,
        pressurant_gas_mass_kg: float | None = None,
        previous_state: TankState | None = None,
        phase_override: str | None = None
    ) -> TankState:

        # if a phase is specified then use it, otherwise use the tank default
        active_model = phase_override if phase_override is not None else self.phase_model

        if active_model == "self_pressurised":
            if total_mass_kg is None:
                raise ValueError("self_pressurised model requires total_mass_kg")
            
            return self._state_from_mass_and_energy_self_pressurised(
                total_mass_kg=total_mass_kg,
                total_internal_energy_j=total_internal_energy_j,
                previous_state=previous_state
            )


        if active_model == "single_phase":
            if total_mass_kg is None:
                raise ValueError("single_phase model requires total_mass_kg")
            
            return self._state_from_mass_and_energy_single_phase(
                total_mass_kg=total_mass_kg,
                total_internal_energy_j=total_internal_energy_j,
                previous_state=previous_state
            )

        if active_model == "pressurised_liquid":
            if liquid_mass_kg is None:
                raise ValueError("pressurised_liquid model requires liquid_mass_kg")
            if pressurant_gas_mass_kg is None:
                raise ValueError("pressurised_liquid model requires pressurant_gas_mass_kg")
            
            return self._state_from_mass_and_energy_pressurised_liquid(
                liquid_mass_kg=liquid_mass_kg,
                pressurant_gas_mass_kg=pressurant_gas_mass_kg,
                total_internal_energy_j=total_internal_energy_j,
                previous_state=previous_state
            )


        raise NotImplementedError(f"Unsupported phase model: {active_model}")


    


@dataclass
class TankState:
    config: TankConfig

    pressure_pa: float | None = None
    temperature_k: float | None = None
    
    total_mass_kg: float | None = None
    total_internal_energy_j: float | None = None

    liquid_mass_kg: float | None = None
    vapour_mass_kg: float | None = None
    pressurant_gas_mass_kg: float | None = None

    ullage_volume_m3: float | None = None
    liquid_volume_m3: float | None = None

    pressurant_gas_internal_energy_j: float | None = None




