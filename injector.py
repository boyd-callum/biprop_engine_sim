

from dataclasses import dataclass
from typing import Literal


# -------------------------------
# Dataclasses for the injectors


@dataclass
class InjectorConfig:
    cd: float
    area_m2: float
