from __future__ import annotations


from dataclasses import dataclass
from typing import Literal
from CoolProp.CoolProp import PropsSI
from constants import *



# -------------------------------
# Dataclasses for the fluids



CoolPropBackend = Literal["HEOS", "REFPROP", "INCOMP", "IF97"]

@dataclass
class Fluid:
    name: str
    coolprop_name: str
    backend: CoolPropBackend = "HEOS"

    cea_name: str | None = None


    @property
    def coolprop_key(self) -> str:
        return f"{self.backend}::{self.coolprop_name}"
        # return self.coolprop_name
    
    
    def props_si(
        self,
        output: str,
        input_1_name: str,
        input_1_value: float,
        input_2_name: str,
        input_2_value: float,
    ) -> float:
        return PropsSI(
            output,
            input_1_name,
            input_1_value,
            input_2_name,
            input_2_value,
            self.coolprop_key,
        )
    

    def get_molar_mass(self) -> float:
        return PropsSI("MOLARMASS", self.coolprop_key)
    

    def get_R(self) -> float:
        # specific gas constant
        return R_UNIVERSAL / self.get_molar_mass()
    

    def get_saturated_vapour_gamma(self, temperature_k: float) -> float:
        """
        Return gamma = cp/cv for saturated vapour at the given temperature.
        """

        cp = PropsSI("Cpmass", "T", temperature_k, "Q", 1.0, self.coolprop_key)
        cv = PropsSI("Cvmass", "T", temperature_k, "Q", 1.0, self.coolprop_key)

        return cp / cv


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

        if (
            saturation_pressure_pa is not None
            and pressure_pa >= saturation_pressure_pa * (1.0 - SATURATION_PRESSURE_TOLERANCE)
        ):
            return self.get_saturated_vapour_gamma(temperature_k)

        try:
            # ratio of specific heats at some given pressure and temperature
            cp = PropsSI("Cpmass", "P", pressure_pa, "T", temperature_k, self.coolprop_key)
            cv = PropsSI("Cvmass", "P", pressure_pa, "T", temperature_k, self.coolprop_key)

            return cp / cv

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
        return PropsSI("Tmin", self.coolprop_key)

    def get_Tmax(self) -> float:
        # highest temperature the fluid model is valid at
        return PropsSI("Tmax", self.coolprop_key)

    def get_Ttriple(self) -> float:
        # returns the triple point temperature
        return PropsSI("Ttriple", self.coolprop_key)
    
    def get_Tcrit(self) -> float:
        # returns critical temperature
        return PropsSI("Tcrit", self.coolprop_key)
    
    def get_Ptriple(self) -> float:
        # returns the triple point pressure
        return PropsSI("ptriple", self.coolprop_key)

    def get_Pcrit(self) -> float:
        # returns critical pressure
        return PropsSI("pcrit", self.coolprop_key)

    def get_saturation_properties_from_temp(self, T: float):
        # returns saturated liquid/vapour properties at a given temp (valid below Tcrit)

        Ttriple = self.get_Ttriple()
        Tcrit = self.get_Tcrit()

        if T > Ttriple and T < Tcrit:
            return dict(
                psat=PropsSI("P", "T", T, "Q", 0.0, self.coolprop_key),
                tsat = T,
                vf=1.0 / PropsSI("D", "T", T, "Q", 0.0, self.coolprop_key),
                vg=1.0 / PropsSI("D", "T", T, "Q", 1.0, self.coolprop_key),
                uf=PropsSI("U", "T", T, "Q", 0.0, self.coolprop_key),
                ug=PropsSI("U", "T", T, "Q", 1.0, self.coolprop_key),
                hf=PropsSI("H", "T", T, "Q", 0.0, self.coolprop_key),
                hg=PropsSI("H", "T", T, "Q", 1.0, self.coolprop_key),
            )
        else:
            raise ValueError(f"{self.name} temperature ({T} K) is not between the Triple point ({Ttriple} K) and Critical Point ({Tcrit} K).")

    def get_saturation_properties_from_pressure(self, P: float):

        Ptriple = self.get_Ptriple()
        Pcrit = self.get_Pcrit()

        if P > Ptriple and P < Pcrit:
            Tsat = PropsSI("T", "P", P, "Q", 0.0, self.coolprop_key)
            sat_props = self.get_saturation_properties_from_temp(Tsat)
            return sat_props

        else:
            raise ValueError(f"{self.name} pressure ({P} Pa) is not between the Triple point ({Ptriple} Pa) and Critical Point ({Pcrit} Pa).")
    
    def get_liquid_enthalpy_from_pressure(self, P: float) -> float:
        """
        returns enthalpy of fluid in the liquid phase, in J/kg
        """
        liquid_enthalpy_j_kg = PropsSI("H", "P", P, "Q", 0.0, self.coolprop_key)

        return liquid_enthalpy_j_kg
    
    def get_vapour_enthalpy_from_pressure(self, P: float) -> float:
        """
        returns enthalpy of fluid in the vapour phase, in J/kg
        """
        vapour_enthalpy_j_kg = PropsSI("H", "P", P, "Q", 1.0, self.coolprop_key)

        return vapour_enthalpy_j_kg
    
    def get_liquid_enthalpy_from_temperature(self, T: float) -> float:
        """
        returns enthalpy of fluid in the liquid phase, in J/kg
        """
        liquid_enthalpy_j_kg = PropsSI("H", "T", T, "Q", 0.0, self.coolprop_key)

        return liquid_enthalpy_j_kg
    
    def get_vapour_enthalpy_from_temperature(self, T: float) -> float:
        """
        returns enthalpy of fluid in the vapour phase, in J/kg
        """
        vapour_enthalpy_j_kg = PropsSI("H", "T", T, "Q", 1.0, self.coolprop_key)

        return vapour_enthalpy_j_kg
    
    def get_fluid_density_from_pressure_temperature(self, P: float, T: float) -> float:
        """
        returns density of a fluid from the pressure and temp, kg/m^3
        """
        density_kg_m3 = PropsSI("D", "T", T, "P", P, self.coolprop_key)

        return density_kg_m3
    
    def get_specific_internal_energy_from_pressure_temperature(self, P: float, T: float) -> float:
        """"
        returns the specific internal energy of a fluid from the pressure and temp, J/kg
        """

        specific_internal_energy_j_kg = PropsSI("U", "P", P, "T", T, self.coolprop_key)

        return specific_internal_energy_j_kg
    
    def get_fluid_enthalpy_from_pressure_temperature(self, P: float, T: float) -> float:
        """
        returns the specific enthalpy of a fluid from the pressure and temp, J/kg
        """
        enthalpy_j_kg = PropsSI("H", "P", P, "T", T, self.coolprop_key)

        return enthalpy_j_kg
    
    def get_fluid_enthalpy_from_density_temperature(self, D: float, T: float) -> float:
        """
        returns the specific enthalpy of a fluid from the density and temp, J/kg
        """
        enthalpy_j_kg = PropsSI("H", "D", D, "T", T, self.coolprop_key)

        return enthalpy_j_kg