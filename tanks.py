from __future__ import annotations


from dataclasses import dataclass
from typing import Literal, Callable

from constants import *
from fluid import Fluid
from helpers import brent_search

# -------------------------------
# Solver settings

# the reconstructed vapour must fill the ullage to within this
SATURATED_VOLUME_TOLERANCE_M3 = 1e-7

# initial half-width of the temperature bracket, as a fraction of the guess
TEMPERATURE_BRACKET_FRACTION = 0.05

# how far the bracket is allowed to widen before giving up
TEMPERATURE_BRACKET_MAX_EXPANSIONS = 12

# how far inside the dome the saturated bracket is held - properties are
# undefined exactly at the triple and critical points
SATURATION_BRACKET_MARGIN_K = 0.01

# convergence of the dew-point temperature solve, in kelvin
DEW_TEMPERATURE_TOLERANCE_K = 1e-9

# convergence of the saturated-tank energy solve, in joules
SATURATED_ENERGY_TOLERANCE_J = 1e-3

# convergence of the pressurised-liquid energy solve, in joules
PRESSURISED_LIQUID_ENERGY_TOLERANCE_J = 1e-2

# convergence of the ullage pressure fixed point, in pascals
PRESSURE_FIXED_POINT_TOLERANCE_PA = 1e-3

# passes of the pressure fixed point before falling back to a bracketed search
PRESSURE_FIXED_POINT_MAX_ITERATIONS = 20

# convergence of that fallback search, as an unfilled tank volume
PRESSURE_RESIDUAL_TOLERANCE_M3 = 1e-9


def solve_temperature_from_guess(
    residual_func: Callable[[float], float],
    temperature_guess_k: float,
    minimum_temperature_k: float,
    maximum_temperature_k: float,
    tolerance: float = 1e-3,
    tank_name: str = "tank",
) -> float:
    """
    Solve for the temperature that zeroes a residual, bracketing outward from a
    guess and never evaluating outside the given valid range.

    Growing the bracket from the previous temperature rather than starting from
    the full valid range matters for two reasons:

      - CoolProp extrapolates silently outside its range instead of raising, so a
        solver given wide fixed bounds can converge on physically meaningless
        states.
      - The residuals here are not globally monotonic. The saturated energy in
        particular turns over near the critical point once the phase split clamps,
        so the two ends of the full range can share a sign even though a root
        exists near the previous temperature. A bracketing solver evaluates its
        endpoints and would refuse the whole solve.

    Both of those are invisible to a plain bisection that only ever probes
    midpoints, which is how the original code got away with them.
    """

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
                tolerance=tolerance,
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


    def get_dew_temperature_k(
            self,
            total_mass_kg: float
    ) -> float:
        """
        The temperature at which this tank's bulk specific volume meets the
        saturated-vapour line - the top of the dome for this fill level.

        Above it the tank is superheated vapour and the two-phase model does not
        apply, so this is the upper bound of any saturated solve. Bounding the
        search here is what keeps it single-rooted: the clamped energy carries a
        second, spurious root in the all-vapour region above this temperature,
        where the state is not on the dome at all.
        """

        bulk_specific_volume_m3_kg = self.tank_volume_m3 / total_mass_kg

        lowest_temperature_k = self.fluid.get_Ttriple() + SATURATION_BRACKET_MARGIN_K
        highest_temperature_k = self.fluid.get_Tcrit() - SATURATION_BRACKET_MARGIN_K

        def vapour_volume_residual(temperature_k: float) -> float:
            saturation_properties = self.fluid.get_saturation_properties_from_temp(temperature_k)
            return saturation_properties['vg'] - bulk_specific_volume_m3_kg

        # vg falls monotonically with temperature, so this has at most one root
        if vapour_volume_residual(highest_temperature_k) > 0.0:
            # saturated vapour is still less dense than the tank right up to the
            # critical point - the dome is available across the whole range
            return highest_temperature_k

        if vapour_volume_residual(lowest_temperature_k) < 0.0:
            raise ValueError(
                f"{self.name} bulk specific volume "
                f"{bulk_specific_volume_m3_kg:.6e} m^3/kg is below the saturated "
                f"vapour line across the whole dome - no two-phase state exists."
            )

        return brent_search(
            residual_func=vapour_volume_residual,
            bounds=[lowest_temperature_k, highest_temperature_k],
            tolerance=0.0,
            x_tolerance=DEW_TEMPERATURE_TOLERANCE_K,
        )


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

        def energy_residual_j(temperature_k: float) -> float:
            return (
                self.calculate_internal_energy(total_mass_kg, temperature_k)
                - total_internal_energy_j
            )

        # Find the saturation temperature that carries the right internal energy at
        # this mass. Only the dome is valid, so that is the bracket.
        #
        # This was a hand-rolled bisection that halved a ~130 K span down to an
        # absolute 1e-3 J tolerance on an energy of order 1e6 J - roughly 40
        # iterations, each of which is a full saturation property lookup. Brent
        # gets the same answer in well under half that.
        #
        # Both bounds are held just inside the dome. Saturation properties are
        # undefined AT the triple point and the critical point, and unlike the old
        # bisection - which only ever evaluated midpoints - a bracketing solver
        # evaluates the endpoints themselves.
        # Bound the search by the dome itself: from just above the triple point up
        # to the dew temperature for this fill level. Inside that range the energy
        # is monotonic in temperature and there is exactly one root.
        #
        # The old bisection used the full triple-to-critical span. That span also
        # contains a second, spurious root above the dew point, where the phase
        # split has clamped to all-vapour and the energy follows ug(T) back down
        # through the target. Bisection happened to converge on the physical root
        # from its starting side; anything that examines the endpoints does not.
        lowest_temperature_k = self.fluid.get_Ttriple() + SATURATION_BRACKET_MARGIN_K
        dew_temperature_k = self.get_dew_temperature_k(total_mass_kg)

        temperature_k = brent_search(
            residual_func=energy_residual_j,
            bounds=[lowest_temperature_k, dew_temperature_k],
            tolerance=SATURATED_ENERGY_TOLERANCE_J,
        )

        correct_temperature_k = temperature_k

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
        temperature_k = solve_temperature_from_guess(
            residual_func=temperature_residual,
            temperature_guess_k=temperature_guess_k,
            minimum_temperature_k=self.fluid.get_Tmin(),
            maximum_temperature_k=self.fluid.get_Tmax(),
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

        # Seed for the pressure fixed point. The temperature solve calls it many
        # times at nearby temperatures, so carrying the last answer forward means
        # it usually starts within a few hundred Pa of the root.
        if previous_state is not None and previous_state.pressure_pa is not None:
            previous_pressure_guess_pa = [previous_state.pressure_pa]
        else:
            previous_pressure_guess_pa = [ATMOSPHERE_PRESSURE_PA * 10.0]


        def pressure_residual_for_temperature(
                pressure_pa: float,
                temperature_k: float
            ) -> float:
            """
            Unfilled tank volume at this pressure and temperature. Zero when the
            liquid and the ullage gas exactly fill the tank.
            """

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


        def solve_pressure_for_temperature(
                temperature_k: float
            ) -> float:

            """
            With a given temperature, find the pressure at which the liquid and the
            ullage gas exactly fill the tank.

            This used to be a bisection over [1e4, 1e8] Pa to a 1e-9 m^3 tolerance -
            around 24 iterations, each costing two real-fluid density evaluations,
            run once per iteration of the temperature solve above it. That nested
            pair was about 80% of the whole simulation's runtime.

            It does not need a search. The volume constraint can be inverted
            directly: the liquid volume fixes the ullage volume, the ullage volume
            and the pressurant mass fix the gas density, and the gas density with
            the temperature gives the pressure. Only the liquid density depends on
            the pressure, and for a liquid that dependence is very weak, so
            iterating on it converges in two or three passes.
            """

            if self.pressurant_fluid is None:
                raise ValueError

            pressure_pa = previous_pressure_guess_pa[0]

            for _ in range(PRESSURE_FIXED_POINT_MAX_ITERATIONS):

                liquid_density_kg_m3 = self.fluid.get_fluid_density_from_pressure_temperature(
                    P=pressure_pa,
                    T=temperature_k
                )

                ullage_volume_m3 = self.tank_volume_m3 - liquid_mass_kg / liquid_density_kg_m3

                if ullage_volume_m3 <= 0.0:
                    # liquid alone over-fills the tank at this pressure - fall back
                    # to the bracketed search, which can cope with it
                    break

                gas_density_kg_m3 = pressurant_gas_mass_kg / ullage_volume_m3

                new_pressure_pa = self.pressurant_fluid.props_si(
                    "P", "D", gas_density_kg_m3, "T", temperature_k
                )

                converged = (
                    abs(new_pressure_pa - pressure_pa)
                    <= PRESSURE_FIXED_POINT_TOLERANCE_PA
                )

                pressure_pa = new_pressure_pa

                if converged:
                    previous_pressure_guess_pa[0] = pressure_pa
                    return pressure_pa

            # Did not settle - fall back to the original bracketed search so a hard
            # case still gets an answer rather than a half-converged guess.
            pressure_pa = brent_search(
                residual_func=lambda p: pressure_residual_for_temperature(p, temperature_k),
                bounds=[1e4, 1e8],
                tolerance=PRESSURE_RESIDUAL_TOLERANCE_M3
            )

            previous_pressure_guess_pa[0] = pressure_pa

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

        temperature_k = brent_search(
            residual_func=temperature_residual,
            bounds=[temp_low_k, temp_high_k],
            tolerance=PRESSURISED_LIQUID_ENERGY_TOLERANCE_J
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




