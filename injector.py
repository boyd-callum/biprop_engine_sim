
from fluid import Fluid
from tanks import TankState

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
            tank_state: TankState,
            downstream_pressure_pa: float
    ) -> float:
        """
        liquid flow through an injector, using the incompressible orifice flow equation
        mdot = cd * A * sqrt(2 * rho * delta_p)
        """


        # get upstream conditions
        upstream_pressure_pa = tank_state.pressure_pa
        upstream_temperature_k = tank_state.temperature_k

        if upstream_pressure_pa is None or upstream_temperature_k is None:
            raise ValueError("Tank state is missing pressure or temperature.")

        liquid_density_kg_m3 = tank_state.config.fluid.get_saturation_properties_from_temp(upstream_temperature_k)['vf']**-1

        # calculate mass flow rate using incompressible orifice flow equation

        delta_p_pa = upstream_pressure_pa - downstream_pressure_pa

        if delta_p_pa <= 0.0:
            raise ValueError(f"Reversed flow in injector. Delta p of {delta_p_pa} Pa")
        
        mdot_kg_s = self.cd * self.area_m2 * math.sqrt(2.0 * liquid_density_kg_m3 * delta_p_pa)

        return mdot_kg_s
    

    def get_gas_mdot_kg_s(
            self,
            tank_state: TankState,
            downstream_pressure_pa: float,
    ) -> float:
        
        """
        vapour-only flow through an injector, using ideal-gas isentropic flow

        choked:
            mdot = cd * A * P0 * sqrt(gamma/R/T) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

        unchoked:
            standard compressible orifice flow
            mdot = 
         
        """

        if tank_state.pressure_pa is None or tank_state.temperature_k is None:
            raise ValueError("Tank state is missing pressure or temperature.")
        
        # get upstream conditions

        upstream_pressure_pa = tank_state.pressure_pa
        upstream_temperature_k = tank_state.temperature_k

        gamma = tank_state.config.fluid.get_gamma_at_PT(upstream_pressure_pa, upstream_temperature_k)

        # determine if flow is choked

        pressure_ratio = downstream_pressure_pa / upstream_pressure_pa
        critical_pressure_ratio = (2/(gamma+1))**(gamma/(gamma-1))

        if pressure_ratio <= critical_pressure_ratio:
            # choked flow
            R = tank_state.config.fluid.get_R()
            mdot_kg_s = self.cd * self.area_m2 * upstream_pressure_pa * math.sqrt(gamma/R/upstream_temperature_k) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
        
        else:
            # unchoked flow
            R = tank_state.config.fluid.get_R()
            mdot_kg_s = self.cd * self.area_m2 * math.sqrt(2*gamma/(gamma-1)*R*upstream_temperature_k) * (pressure_ratio**(2/gamma) - pressure_ratio**((gamma+1)/gamma))**0.5


        return mdot_kg_s