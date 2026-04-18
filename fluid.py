

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
    

    def get_saturation_properties(self, T: float):
        # returns saturated liquid/vapour properties at a given temp (valid below Tcrit)

        return dict(
            psat=PropsSI("P", "T", T, "Q", 0.0, self.coolprop_key),
            vf=1.0 / PropsSI("D", "T", T, "Q", 0.0, self.coolprop_key),
            vg=1.0 / PropsSI("D", "T", T, "Q", 1.0, self.coolprop_key),
            uf=PropsSI("U", "T", T, "Q", 0.0, self.coolprop_key),
            ug=PropsSI("U", "T", T, "Q", 1.0, self.coolprop_key),
            hf=PropsSI("H", "T", T, "Q", 0.0, self.coolprop_key),
            hg=PropsSI("H", "T", T, "Q", 1.0, self.coolprop_key),
        )
    
