

from dataclasses import dataclass
from rocketcea.cea_obj_w_units import CEA_Obj
from typing import Literal

from constants import ATMOSPHERE_PRESSURE_PA, STANDARD_GRAVITY_M_S2, AMBIENT_TEMPERATURE_K

# -------------------------------
# Dataclasses for the engine

@dataclass
class EngineGeometry:
    nozzle_throat_area_m2: float
    expansion_ratio: float

@dataclass
class EngineConfig:
    geometry: EngineGeometry
    state: EngineState | None = None
    ox_cea_name: str | None = None
    fuel_cea_name: str | None = None
    cstar_efficiency: float = 1.0
    cf_efficiency: float = 1.0
    

    def __post_init__(self) -> None:
        self.state = EngineState(
            model=self,
            # chamber_pressure_pa=ATMOSPHERE_PRESSURE_PA,
            chamber_pressure_pa=30 * 1e5, # 30 bar
            chamber_temperature_k=AMBIENT_TEMPERATURE_K,
        )


    def make_cea_obj(self) -> CEA_Obj:
        
        if self.ox_cea_name is None or self.fuel_cea_name is None:
            raise ValueError("ox_cea_name or fuel_cea_name not defined.")

        return CEA_Obj(
            oxName=self.ox_cea_name,
            fuelName=self.fuel_cea_name,
            pressure_units="Pa",
            cstar_units="m/s",
            temperature_units="K",
            sonic_velocity_units="m/s",
            enthalpy_units="J/kg",
            density_units="kg/m^3",
            specific_heat_units="J/kg-K"
        )




    def calculate_state(
        self,
        ox_mdot_kg_s: float,
        fuel_mdot_kg_s: float,
        ambient_pressure_pa: float,
        previous_state: EngineState
   ) -> EngineState:
        

        min_mdot_kg_s = 1e-9

        if ox_mdot_kg_s <= min_mdot_kg_s or fuel_mdot_kg_s <= min_mdot_kg_s:
            return EngineState(
                model=self,
                ox_mdot_kg_s=ox_mdot_kg_s,
                fuel_mdot_kg_s=fuel_mdot_kg_s,
                total_mdot_kg_s=ox_mdot_kg_s + fuel_mdot_kg_s,
                mixture_ratio=None,
                chamber_pressure_pa=ambient_pressure_pa,
                chamber_temperature_k=AMBIENT_TEMPERATURE_K,
                cstar_ideal_m_s=None,
                cstar_delivered_m_s=None,
                cf_ideal=None,
                cf_delivered=None,
                thrust_n=0.0,
                isp_s=0.0,
            )


        previous_chamber_pressure_pa = previous_state.chamber_pressure_pa
        
        mixture_ratio = ox_mdot_kg_s / fuel_mdot_kg_s
        total_mdot_kg_s = ox_mdot_kg_s + fuel_mdot_kg_s

        cea = self.make_cea_obj()

        # inital chamber pressure guess
        chamber_pressure_pa = previous_chamber_pressure_pa or ATMOSPHERE_PRESSURE_PA

        pc_tolerance = 1e-4

        if self.state is not None:
            cstar_ideal_m_s = self.state.cstar_ideal_m_s
            cstar_delivered_m_s = self.state.cstar_delivered_m_s
        else:
            cstar_ideal_m_s = None
            cstar_delivered_m_s = None

        # theres a circular relationship between cstar and chmaber pressure, so we do a couple iterations to find a more accurate guess
        for _ in range(3):
            cstar_ideal_m_s = cea.get_Cstar(Pc=chamber_pressure_pa, MR=mixture_ratio)
            cstar_delivered_m_s = self.cstar_efficiency * cstar_ideal_m_s

            new_chamber_pressure_pa = total_mdot_kg_s * cstar_delivered_m_s / self.geometry.nozzle_throat_area_m2
            

            if abs(new_chamber_pressure_pa - chamber_pressure_pa) / chamber_pressure_pa < pc_tolerance:
                break
            
            chamber_pressure_pa = new_chamber_pressure_pa

        chamber_temperature_k = cea.get_Tcomb(
            Pc=chamber_pressure_pa,
            MR=mixture_ratio
        )

        # finding thrust coeff
        cf_result = cea.get_PambCf(
            Pamb=ambient_pressure_pa,
            Pc=chamber_pressure_pa,
            MR=mixture_ratio,
            eps=self.geometry.expansion_ratio
        )

        cf_ideal = cf_result[1]
        cf_delivered = self.cf_efficiency * cf_ideal

        thrust_n = cf_delivered * chamber_pressure_pa * self.geometry.nozzle_throat_area_m2

        isp_s = thrust_n / (total_mdot_kg_s * STANDARD_GRAVITY_M_S2)
        

        return EngineState(
            model=self,
            ox_mdot_kg_s=ox_mdot_kg_s,
            fuel_mdot_kg_s=fuel_mdot_kg_s,
            total_mdot_kg_s=total_mdot_kg_s,
            mixture_ratio=mixture_ratio,
            chamber_pressure_pa=chamber_pressure_pa,
            chamber_temperature_k=chamber_temperature_k,
            cstar_ideal_m_s=cstar_ideal_m_s,
            cstar_delivered_m_s=cstar_delivered_m_s,
            cf_ideal=cf_ideal,
            cf_delivered=cf_delivered,
            thrust_n=thrust_n,
            isp_s=isp_s,
        )
        




@dataclass
class EngineState:
    model: EngineConfig

    chamber_pressure_pa: float
    chamber_temperature_k: float

    ox_mdot_kg_s: float | None = None
    fuel_mdot_kg_s: float | None = None
    total_mdot_kg_s: float | None = None

    mixture_ratio: float | None = None

    cstar_ideal_m_s: float | None = None
    cstar_delivered_m_s: float | None = None

    cf_ideal: float | None = None
    cf_delivered: float | None = None

    thrust_n: float | None = None
    isp_s: float | None = None
