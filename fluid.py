

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
    

    def get_gamma_at_PT(self, pressure_pa: float, temperature_k: float) -> float:
        """"
        Return gamma = cp/cv at given temp and pressure

        For nitrous, coolprop can fail if the state lies exactly on or close to the saturation curve, so this function nudges the state slightly into the vapour region before retrying
        """
        try:
            # ratio of specific heats at some given pressure and temperature
            cp = PropsSI("Cpmass", "P", pressure_pa, "T", temperature_k, self.coolprop_key)
            cv = PropsSI("Cvmass", "P", pressure_pa, "T", temperature_k, self.coolprop_key)
            
            gamma = cp / cv

            return gamma
        except Exception as original_error:

            try:
                
                psat = self.get_saturation_properties_from_temp(temperature_k)["psat"]

                # if we are on/above saturation, nudge slightly into vapour region
                nudged_pressure_pa = min(pressure_pa, 0.999 * psat)

                cp = PropsSI("Cpmass", "P", nudged_pressure_pa, "T", temperature_k, self.coolprop_key)
                cv = PropsSI("Cvmass", "P", nudged_pressure_pa, "T", temperature_k, self.coolprop_key)

                gamma = cp / cv

                return gamma

            except Exception as fallback_error:
                raise RuntimeError(
                    "get_gamma_at_PT failed. Inputs: \n"
                    f"Pressure: {pressure_pa:.2f} Pa\n"
                    f"Temperature: {temperature_k:.2f} K\n"
                    f"Original error: {original_error}\n"
                    f"Fallback error: {fallback_error}"
                    )

    


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