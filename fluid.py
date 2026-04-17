

from dataclasses import dataclass
from typing import Literal
from CoolProp.CoolProp import PropsSI



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