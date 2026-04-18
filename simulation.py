
from dataclasses import dataclass
from typing import Literal


from tanks import TankState
from engine import EngineState
from cases import SimCase, n2o_blowdown_case

# -------------------------------
# Dataclasses for the simulation running

 


@dataclass
class SimPoint:
    """A point in the simulation, with a timestamp and states for the engine and fluid sources."""
    time_s: float

    engine: EngineState
    tanks: dict[str, TankState]
    


def blowdown_simulate(
        case: SimCase,
        record: bool = True
):
    return None