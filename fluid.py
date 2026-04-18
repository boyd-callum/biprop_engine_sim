

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
    

    def get_gamma_at_PT(self, P: float, T: float) -> float:
        # ratio of specific heats at some given pressure and temperature
        cp = PropsSI("Cpmass", "P", P, "T", T, self.coolprop_key)
        cv = PropsSI("Cvmass", "P", P, "T", T, self.coolprop_key)

        return cp / cv
    
    def get_Ttriple(self) -> float:
        # returns the triple point temperature
        return PropsSI("Ttriple", self.coolprop_key)
    
    def get_Tcrit(self) -> float:
        # returns critical temperature
        return PropsSI("Tcrit", self.coolprop_key)
    
    def get_Ptriple(self) -> float:
        # returns the triple point pressure
        return PropsSI("Ptriple", self.coolprop_key)

    def get_Pcrit(self) -> float:
        # returns critical pressure
        return PropsSI("Pcrit", self.coolprop_key)

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