
from dataclasses import dataclass
from typing import Literal


from tanks import TankState
from engine import EngineState

# -------------------------------
# Dataclasses for the simulation running

 


@dataclass
class SimPoint:
    """A point in the simulation, with a timestamp and states for the engine and fluid sources."""
    time_s: float

    engine: EngineState
    tanks: dict[str, TankState]
    