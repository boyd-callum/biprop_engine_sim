from __future__ import annotations


from dataclasses import dataclass, field
from typing import Literal
import CoolProp as CP
from CoolProp.CoolProp import PropsSI, AbstractState
from constants import *



# -------------------------------
# Dataclasses for the fluids



CoolPropBackend = Literal["HEOS", "REFPROP", "INCOMP", "IF97"]


# Input-pair codes for the low-level interface, keyed by the (name_1, name_2)
# strings the PropsSI-style calls use. QT and PQ take their arguments in the
# opposite order to the string form, which _input_pair_arguments handles.
_INPUT_PAIRS = {
    ("D", "T"): CP.DmassT_INPUTS,
    ("T", "D"): CP.DmassT_INPUTS,
    ("P", "T"): CP.PT_INPUTS,
    ("T", "P"): CP.PT_INPUTS,
    ("P", "Q"): CP.PQ_INPUTS,
    ("Q", "P"): CP.PQ_INPUTS,
    ("T", "Q"): CP.QT_INPUTS,
    ("Q", "T"): CP.QT_INPUTS,
    ("P", "S"): CP.PSmass_INPUTS,
    ("S", "P"): CP.PSmass_INPUTS,
    ("P", "H"): CP.HmassP_INPUTS,
    ("H", "P"): CP.HmassP_INPUTS,
}

# The order each input pair expects its two values in
_INPUT_PAIR_ORDER = {
    CP.DmassT_INPUTS: ("D", "T"),
    CP.PT_INPUTS: ("P", "T"),
    CP.PQ_INPUTS: ("P", "Q"),
    CP.QT_INPUTS: ("Q", "T"),
    CP.PSmass_INPUTS: ("P", "S"),
    CP.HmassP_INPUTS: ("H", "P"),
}

# Output accessors on the low-level state object
_OUTPUTS = {
    "D": "rhomass",
    "H": "hmass",
    "P": "p",
    "Q": "Q",
    "S": "smass",
    "T": "T",
    "U": "umass",
    "Cpmass": "cpmass",
    "Cvmass": "cvmass",
}


@dataclass
class Fluid:
    name: str
    coolprop_name: str
    backend: CoolPropBackend = "HEOS"

    cea_name: str | None = None

    # Lazily built low-level state object and the fixed-point constants.
    # Neither is part of the fluid definition, so both stay out of comparison.
    _state: object | None = field(default=None, init=False, repr=False, compare=False)
    _constants: dict = field(default_factory=dict, init=False, repr=False, compare=False)
    _key: str | None = field(default=None, init=False, repr=False, compare=False)


    @property
    def coolprop_key(self) -> str:
        # built once - this was previously an f-string rebuilt on every property
        # call, which ran to hundreds of thousands of times per simulation
        if self._key is None:
            self._key = f"{self.backend}::{self.coolprop_name}"
        return self._key


    def get_state(self):
        """
        The low-level CoolProp state object for this fluid.

        PropsSI re-parses the fluid string and rebuilds the backend on every call,
        which costs about 95 us. Reusing one AbstractState and calling update() on
        it costs about 8 us for the same query - a bit under 12x.

        The object is stateful: every accessor here updates it and reads the result
        immediately, so no caller can observe a half-updated state.
        """

        if self._state is None:
            self._state = AbstractState(self.backend, self.coolprop_name)

        return self._state


    def _update(self, name_1: str, value_1: float, name_2: str, value_2: float):
        """
        Update the low-level state from an input pair, returning it for reading.
        """

        pair = _INPUT_PAIRS.get((name_1, name_2))

        if pair is None:
            raise KeyError(
                f"No low-level input pair for ({name_1}, {name_2}). "
                f"Add it to _INPUT_PAIRS or use PropsSI directly."
            )

        values = {name_1: value_1, name_2: value_2}
        first_name, second_name = _INPUT_PAIR_ORDER[pair]

        state = self.get_state()
        state.update(pair, values[first_name], values[second_name])

        return state


    def props_si(
        self,
        output: str,
        input_1_name: str,
        input_1_value: float,
        input_2_name: str,
        input_2_value: float,
    ) -> float:

        accessor = _OUTPUTS.get(output)

        if accessor is None or (input_1_name, input_2_name) not in _INPUT_PAIRS:
            # anything not covered by the low-level interface falls back to the
            # original call, so adding a new query can never silently break
            return PropsSI(
                output,
                input_1_name,
                input_1_value,
                input_2_name,
                input_2_value,
                self.coolprop_key,
            )

        state = self._update(input_1_name, input_1_value, input_2_name, input_2_value)

        return getattr(state, accessor)()


    def _get_constant(self, name: str) -> float:
        """
        Fixed points of the fluid. These never change, but were being looked up
        through PropsSI tens of thousands of times per run.
        """

        if name not in self._constants:
            self._constants[name] = PropsSI(name, self.coolprop_key)

        return self._constants[name]


    def get_molar_mass(self) -> float:
        return self._get_constant("MOLARMASS")


    def get_R(self) -> float:
        # specific gas constant
        return R_UNIVERSAL / self.get_molar_mass()
    

    def get_saturated_vapour_gamma(self, temperature_k: float) -> float:
        """
        Return gamma = cp/cv for saturated vapour at the given temperature.
        """

        state = self._update("T", temperature_k, "Q", 1.0)

        return state.cpmass() / state.cvmass()


    def get_gamma_at_PT(self, pressure_pa: float, temperature_k: float) -> float:
        """"
        Return gamma = cp/cv at given temp and pressure

        A saturated tank sits exactly on the saturation curve, where pressure and
        temperature are not independent and CoolProp cannot resolve a state from
        the pair. That is a normal operating point for a blowdown tank, not an
        error, so it is detected up front and evaluated as saturated vapour -
        which is the correct property for vapour leaving the ullage anyway.
        """

        # is this state on (or above) the saturation line for this temperature?
        try:
            saturation_pressure_pa = self.get_saturation_properties_from_temp(temperature_k)["psat"]
        except ValueError:
            saturation_pressure_pa = None

        # Only a state sitting ON the saturation line is ambiguous. A pressure
        # merely above the saturation pressure is compressed liquid and a pressure
        # below it is superheated vapour - both resolve fine from P and T, and
        # treating either as saturated vapour would return the wrong gamma
        # entirely (ethanol at 50 bar, 293 K is liquid, not vapour).
        if (
            saturation_pressure_pa is not None
            and abs(pressure_pa - saturation_pressure_pa)
            <= SATURATION_PRESSURE_TOLERANCE * saturation_pressure_pa
        ):
            return self.get_saturated_vapour_gamma(temperature_k)

        try:
            # ratio of specific heats at some given pressure and temperature
            state = self._update("P", pressure_pa, "T", temperature_k)

            return state.cpmass() / state.cvmass()

        except Exception as original_error:

            if saturation_pressure_pa is None:
                raise RuntimeError(
                    "get_gamma_at_PT failed and the state is not saturated. Inputs: \n"
                    f"Pressure: {pressure_pa:.2f} Pa\n"
                    f"Temperature: {temperature_k:.2f} K\n"
                    f"Original error: {original_error}"
                ) from original_error

            return self.get_saturated_vapour_gamma(temperature_k)

    


    def get_Tmin(self) -> float:
        # lowest temperature the fluid model is valid at
        return self._get_constant("Tmin")

    def get_Tmax(self) -> float:
        # highest temperature the fluid model is valid at
        return self._get_constant("Tmax")

    def get_Ttriple(self) -> float:
        # returns the triple point temperature
        return self._get_constant("Ttriple")

    def get_Tcrit(self) -> float:
        # returns critical temperature
        return self._get_constant("Tcrit")

    def get_Ptriple(self) -> float:
        # returns the triple point pressure
        return self._get_constant("ptriple")

    def get_Pcrit(self) -> float:
        # returns critical pressure
        return self._get_constant("pcrit")

    def get_saturation_properties_from_temp(self, T: float):
        """
        Saturated liquid/vapour properties at a given temperature (below Tcrit).

        This is the hottest property call in the simulation. It used to make nine
        PropsSI calls - seven properties plus the triple and critical points, which
        are constants - and is now two updates of the low-level state, one per
        phase, reading four properties from each.
        """

        if T > self.get_Ttriple() and T < self.get_Tcrit():

            state = self._update("T", T, "Q", 0.0)
            psat = state.p()
            vf = 1.0 / state.rhomass()
            uf = state.umass()
            hf = state.hmass()

            state = self._update("T", T, "Q", 1.0)

            return dict(
                psat=psat,
                tsat=T,
                vf=vf,
                vg=1.0 / state.rhomass(),
                uf=uf,
                ug=state.umass(),
                hf=hf,
                hg=state.hmass(),
            )
        else:
            raise ValueError(f"{self.name} temperature ({T} K) is not between the Triple point ({self.get_Ttriple()} K) and Critical Point ({self.get_Tcrit()} K).")

    def get_saturation_properties_from_pressure(self, P: float):

        Ptriple = self.get_Ptriple()
        Pcrit = self.get_Pcrit()

        if P > Ptriple and P < Pcrit:
            Tsat = self._update("P", P, "Q", 0.0).T()
            sat_props = self.get_saturation_properties_from_temp(Tsat)
            return sat_props

        else:
            raise ValueError(f"{self.name} pressure ({P} Pa) is not between the Triple point ({Ptriple} Pa) and Critical Point ({Pcrit} Pa).")
    
    def get_liquid_enthalpy_from_pressure(self, P: float) -> float:
        """
        returns enthalpy of fluid in the liquid phase, in J/kg
        """
        liquid_enthalpy_j_kg = self._update("P", P, "Q", 0.0).hmass()

        return liquid_enthalpy_j_kg
    
    def get_vapour_enthalpy_from_pressure(self, P: float) -> float:
        """
        returns enthalpy of fluid in the vapour phase, in J/kg
        """
        vapour_enthalpy_j_kg = self._update("P", P, "Q", 1.0).hmass()

        return vapour_enthalpy_j_kg
    
    def get_liquid_enthalpy_from_temperature(self, T: float) -> float:
        """
        returns enthalpy of fluid in the liquid phase, in J/kg
        """
        liquid_enthalpy_j_kg = self._update("T", T, "Q", 0.0).hmass()

        return liquid_enthalpy_j_kg
    
    def get_vapour_enthalpy_from_temperature(self, T: float) -> float:
        """
        returns enthalpy of fluid in the vapour phase, in J/kg
        """
        vapour_enthalpy_j_kg = self._update("T", T, "Q", 1.0).hmass()

        return vapour_enthalpy_j_kg
    
    def get_fluid_density_from_pressure_temperature(self, P: float, T: float) -> float:
        """
        returns density of a fluid from the pressure and temp, kg/m^3
        """
        density_kg_m3 = self._update("P", P, "T", T).rhomass()

        return density_kg_m3
    
    def get_specific_internal_energy_from_pressure_temperature(self, P: float, T: float) -> float:
        """"
        returns the specific internal energy of a fluid from the pressure and temp, J/kg
        """

        specific_internal_energy_j_kg = self._update("P", P, "T", T).umass()

        return specific_internal_energy_j_kg
    
    def get_fluid_enthalpy_from_pressure_temperature(self, P: float, T: float) -> float:
        """
        returns the specific enthalpy of a fluid from the pressure and temp, J/kg
        """
        enthalpy_j_kg = self._update("P", P, "T", T).hmass()

        return enthalpy_j_kg
    
    def get_fluid_enthalpy_from_density_temperature(self, D: float, T: float) -> float:
        """
        returns the specific enthalpy of a fluid from the density and temp, J/kg
        """
        enthalpy_j_kg = self._update("D", D, "T", T).hmass()

        return enthalpy_j_kg