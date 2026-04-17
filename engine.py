

from dataclasses import dataclass
from typing import Literal


# -------------------------------
# Dataclasses for the engine

@dataclass
class EngineGeometry:
    nozzle_throat_area_m2: float
    expansion_ratio: float

@dataclass
class EngineConfig:
    geometry: EngineGeometry
    cstar_efficiency: float = 1.0
    cf_efficiency: float = 1.0

@dataclass
class EngineState:
    model: EngineConfig
    chamber_pressure_pa: float | None = None
    thrust_n: float | None = None
    isp_s: float | None = None
    mixture_ratio: float | None = None
