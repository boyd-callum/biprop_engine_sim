
from fluid import Fluid
from tanks import TankState

from dataclasses import dataclass
from typing import Literal
import math

# -------------------------------
# Dataclasses for the injectors

InjectorRole = Literal["fuel", "oxidiser", "generic"]

@dataclass
class InjectorConfig:
    cd: float
    area_m2: float
    k: float = 2.0
    role: InjectorRole = "generic"

    # flow modeling using dyer method (formally known as the NHNE model)
    # https://wikis.mit.edu/confluence/display/RocketTeam/Modeling


    # weighted average of Homogenous Equilibrium Model (HEM) and Single Phase Incompressible (SPI) models

    def get_SPI_mdot_kg_s(
            self,
            tank_state: TankState,
            downstream_pressure_pa: float
    ) -> float:

        """
        mdot_spi = Cd * A_inj * sqrt(2*rho*delta_P)
        """


        temperature_k = tank_state.temperature_k
        upstream_pressure_pa = tank_state.pressure_pa

        if temperature_k is None or upstream_pressure_pa is None:
            raise ValueError("Tank state is missing pressure or temperature.")

        saturation_properties = tank_state.config.fluid.get_saturation_properties_from_temp(temperature_k)

        density_kg_m3 = 1.0 / saturation_properties['vf']

        delta_p_pa = upstream_pressure_pa - downstream_pressure_pa
        if delta_p_pa <= 0.0:
            return 0.0

        mdot_kg_s = self.cd * self.area_m2 * math.sqrt(2 * density_kg_m3 * delta_p_pa)

        return mdot_kg_s
    

    def get_HEM_mdot_kg_s(
            self,
            tank_state: TankState,
            downstream_pressure_pa: float
    ) -> float:
        
        """"
        mdot_hem = Cd * A_inj * rho_2 * sqrt(2*(h1-h2))

        where:
            - h1 is enthalpy at orifice inlet
            - h2 is enthalpy at orifice outlet
        
        Only valid for liquid flow

        """

        upstream_pressure_pa = tank_state.pressure_pa
        upstream_temperature_k = tank_state.temperature_k

        if upstream_pressure_pa is None or upstream_temperature_k is None:
            raise ValueError("Tank state is missing pressure or temperature")
        
        fluid = tank_state.config.fluid


        # clamp reversed flow
        if downstream_pressure_pa >= upstream_pressure_pa:
            return 0.0

        # HEM is based in liquid inlet state for flashing flow

        # find inlet properties 
        inlet_entropy_j_kg_k = fluid.props_si("S", "P", upstream_pressure_pa, "Q", 0.0)
        inlet_enthalpy_j_kg = fluid.props_si("H", "P", upstream_pressure_pa, "Q", 0.0)

        # find outlet properties
        outlet_density_kg_m3 = fluid.props_si("D", "P", downstream_pressure_pa, "S", inlet_entropy_j_kg_k)
        outlet_enthalpy_j_kg = fluid.props_si("H", "P", downstream_pressure_pa, "S", inlet_entropy_j_kg_k)


        delta_enthalpy = inlet_enthalpy_j_kg - outlet_enthalpy_j_kg

        if delta_enthalpy < 0.0:
            raise ValueError(
                f"Negative HEM enthalpy drop: "
                f"h1={inlet_enthalpy_j_kg} j/kg,    h2={outlet_enthalpy_j_kg} j/kg, "
                f"P1={upstream_pressure_pa}, Pa   P2={downstream_pressure_pa} Pa"
            )

        mdot_kg_s = self.cd * self.area_m2 * outlet_density_kg_m3 * math.sqrt(2.0 * delta_enthalpy)
        
        return mdot_kg_s


    def get_dyer_mdot_kg_s(
            self,
            tank_state: TankState,
            downstream_pressure_pa: float
    ) -> float:
        """"
        dyer is a weighted average of HEM and SPI
        """
        mdot_HEM_kg_s = self.get_HEM_mdot_kg_s(
            tank_state = tank_state,
            downstream_pressure_pa = downstream_pressure_pa
        )

        mdot_SPI_kg_s = self.get_SPI_mdot_kg_s(
            tank_state = tank_state,
            downstream_pressure_pa = downstream_pressure_pa
        )

        k = self.k

        # formula from https://wikis.mit.edu/confluence/display/RocketTeam/Modeling
        mdot_dyer_ks_s = (k * mdot_SPI_kg_s + mdot_HEM_kg_s) / (1.0 + k)

        return mdot_dyer_ks_s


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

        liquid_density_kg_m3 = tank_state.config.fluid.get_fluid_density_from_pressure_temperature(
            P=upstream_pressure_pa,
            T=upstream_temperature_k
        )

        # calculate mass flow rate using incompressible orifice flow equation

        delta_p_pa = upstream_pressure_pa - downstream_pressure_pa

        if delta_p_pa <= 0.0:
            # raise ValueError(f"Reversed flow in injector. Delta p of {delta_p_pa} Pa")
            return 0.0
        
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

        if upstream_pressure_pa - downstream_pressure_pa <= 0.0:
            return 0.0

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
            mdot_kg_s = self.cd * self.area_m2 * upstream_pressure_pa * math.sqrt(2*gamma / ((gamma-1) * R * upstream_temperature_k)) * (pressure_ratio**(2/gamma) - pressure_ratio**((gamma+1)/gamma))**0.5


        return mdot_kg_s