

from dataclasses import dataclass
from typing import Literal
import math

# -------------------------------
# Dataclasses for the injectors


@dataclass
class InjectorConfig:
    cd: float
    area_m2: float


    def get_liquid_mdot_kg_s(
            self,
            upstream_pressure_pa: float,
            downstream_pressure_pa: float, 
            liquid_density_kg_m3: float
    ) -> float:

        delta_p_pa = upstream_pressure_pa - downstream_pressure_pa

        if delta_p_pa <= 0.0:
            raise ValueError(f"Reversed flow in injector. Delta p of {delta_p_pa} Pa")
        
        mdot_kg_s = self.cd * self.area_m2 * math.sqrt(2.0 * liquid_density_kg_m3 * delta_p_pa)

        return mdot_kg_s