from dataclasses import dataclass


@dataclass
class RSEComponent:
    name: str

    length_mm: float
    diameter_mm: float

    dry_mass_g: float

    @property
    def cg_mm(self) -> float:
        """assume CG is center of component"""
        return self.length_mm / 2




