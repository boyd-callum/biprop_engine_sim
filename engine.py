from __future__ import annotations


from dataclasses import dataclass, field
from rocketcea.cea_obj_w_units import CEA_Obj
from typing import Literal, Callable

from constants import ATMOSPHERE_PRESSURE_PA, STANDARD_GRAVITY_M_S2, AMBIENT_TEMPERATURE_K
from helpers import brent_search

# -------------------------------
# Solver settings

# below this, a propellant path is treated as not flowing and the engine as off
MINIMUM_MDOT_KG_S = 1e-9

# passes of the cstar refresh loop in the chamber pressure solve
CHAMBER_PRESSURE_OUTER_ITERATIONS = 4

# convergence of the outer cstar loop, in chamber pressure
CHAMBER_PRESSURE_TOLERANCE_PA = 500.0

# trial pressures within this are treated as the same point by the supply cache
SUPPLY_CACHE_RESOLUTION_PA = 10.0

# convergence of the inner pressure balance, as a mass flow imbalance
CHAMBER_PRESSURE_RESIDUAL_TOLERANCE_KG_S = 1e-5

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

    # RocketCEA objects, keyed by propellant pair. Not part of the engine
    # definition, so excluded from comparison and repr.
    _cea_cache: dict = field(default_factory=dict, init=False, repr=False, compare=False)


    def __post_init__(self) -> None:
        self.reset_state()


    def reset_state(self) -> None:
        """
        Return the engine to its pre-run state.

        Called at the start of every simulation so a run cannot inherit the
        previous run's final chamber pressure - that made results depend on how
        many times the case had already been simulated, which quietly skewed the
        first point of every sweep.
        """

        self.state = EngineState(
            model=self,
            chamber_pressure_pa=ATMOSPHERE_PRESSURE_PA,
            chamber_temperature_k=AMBIENT_TEMPERATURE_K,
        )


    def make_cea_obj(self) -> CEA_Obj:

        if self.ox_cea_name is None or self.fuel_cea_name is None:
            raise ValueError("ox_cea_name or fuel_cea_name not defined.")

        # RocketCEA caches results on the object, keyed by (Pc, MR). Rebuilding it
        # every call throws that cache away, so hold one per propellant pair.
        cache_key = (self.ox_cea_name, self.fuel_cea_name)

        cached = self._cea_cache.get(cache_key)
        if cached is not None:
            return cached

        cea = CEA_Obj(
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

        self._cea_cache[cache_key] = cea

        return cea


    def solve_chamber_pressure_pa(
        self,
        ox_mdot_at_pressure: Callable[[float], float],
        fuel_mdot_at_pressure: Callable[[float], float],
        max_upstream_pressure_pa: float,
        previous_chamber_pressure_pa: float | None = None,
    ) -> float:
        """
        Solve the coupled injector/nozzle pressure balance for chamber pressure.

        The chamber is quasi-steady, so the propellant mass flow the injectors can
        deliver at a given chamber pressure must equal the mass flow the throat can
        pass at that same pressure:

            mdot_ox(Pc) + mdot_fuel(Pc) = Pc * At / cstar(Pc, MR(Pc))

        Injector flow falls with Pc and throat flow rises with it, so the residual
        is monotonic and has exactly one root.

        Solving this implicitly is what keeps the simulation stable. Computing Pc
        explicitly from the previous timestep's mass flow puts a one-step lag in a
        feedback loop whose gain is roughly (cstar / At) * (mdot / 2*dP). That is
        well under one early in a burn, but it crosses one as the tanks drain and
        dP shrinks, at which point the explicit scheme breaks into a period-2
        oscillation between high-flow/high-Pc and low-flow/low-Pc states.
        """

        cea = self.make_cea_obj()
        throat_area_m2 = self.geometry.nozzle_throat_area_m2

        # Absolute bounds: the root cannot be below ambient, and cannot reach the
        # feed pressure because the injectors deliver nothing there.
        floor_pressure_pa = ATMOSPHERE_PRESSURE_PA
        ceiling_pressure_pa = max(max_upstream_pressure_pa, floor_pressure_pa * 1.001)

        chamber_pressure_pa = previous_chamber_pressure_pa or ATMOSPHERE_PRESSURE_PA
        chamber_pressure_pa = min(max(chamber_pressure_pa, floor_pressure_pa), ceiling_pressure_pa)

        # The injector supply depends only on the tank states and the trial
        # pressure, not on cstar, so it is the same across every pass of the outer
        # loop. Each evaluation is a handful of CoolProp calls, so cache them.
        supply_cache: dict[int, tuple[float, float]] = {}

        def get_supply_mdot_kg_s(trial_pressure_pa: float) -> tuple[float, float]:

            cache_key = int(trial_pressure_pa / SUPPLY_CACHE_RESOLUTION_PA)

            cached = supply_cache.get(cache_key)
            if cached is not None:
                return cached

            supply = (
                ox_mdot_at_pressure(trial_pressure_pa),
                fuel_mdot_at_pressure(trial_pressure_pa),
            )

            supply_cache[cache_key] = supply

            return supply

        def make_residual(cstar_delivered_m_s: float):

            def residual_kg_s(trial_pressure_pa: float) -> float:
                """
                Injector supply minus what the throat can pass, at this pressure.
                Supply falls with pressure and throat flow rises with it, so this
                is monotonic and has exactly one root.
                """

                ox_mdot_kg_s, fuel_mdot_kg_s = get_supply_mdot_kg_s(trial_pressure_pa)

                nozzle_mdot_kg_s = trial_pressure_pa * throat_area_m2 / cstar_delivered_m_s

                return ox_mdot_kg_s + fuel_mdot_kg_s - nozzle_mdot_kg_s

            return residual_kg_s

        # cstar is only a weak function of chamber pressure, so hold it fixed while
        # solving the pressure balance and then refresh it. Two or three passes is
        # enough, which keeps CEA calls per timestep in single figures rather than
        # one per root-solve iteration.
        for _ in range(CHAMBER_PRESSURE_OUTER_ITERATIONS):

            ox_mdot_kg_s, fuel_mdot_kg_s = get_supply_mdot_kg_s(chamber_pressure_pa)

            if ox_mdot_kg_s <= MINIMUM_MDOT_KG_S or fuel_mdot_kg_s <= MINIMUM_MDOT_KG_S:
                # no combustion - the engine is not running
                return ATMOSPHERE_PRESSURE_PA

            cstar_delivered_m_s = self.cstar_efficiency * cea.get_Cstar(
                Pc=chamber_pressure_pa,
                MR=ox_mdot_kg_s / fuel_mdot_kg_s
            )

            residual_kg_s = make_residual(cstar_delivered_m_s)

            # Chamber pressure moves by well under a percent per timestep, so start
            # from a narrow bracket around the previous value and widen only if it
            # does not straddle the root. Brent then converges in a few evaluations.
            low_pressure_pa = max(chamber_pressure_pa * 0.8, floor_pressure_pa)
            high_pressure_pa = min(chamber_pressure_pa * 1.25, ceiling_pressure_pa)

            low_residual_kg_s = residual_kg_s(low_pressure_pa)
            high_residual_kg_s = residual_kg_s(high_pressure_pa)

            if low_residual_kg_s < 0.0:
                low_pressure_pa = floor_pressure_pa
                low_residual_kg_s = residual_kg_s(low_pressure_pa)

                if low_residual_kg_s <= 0.0:
                    # not even ambient can be sustained - the engine is not running
                    return ATMOSPHERE_PRESSURE_PA

            if high_residual_kg_s > 0.0:
                high_pressure_pa = ceiling_pressure_pa
                high_residual_kg_s = residual_kg_s(high_pressure_pa)

                if high_residual_kg_s >= 0.0:
                    # would need a chamber pressure at or above the feed pressure
                    return ceiling_pressure_pa

            new_chamber_pressure_pa = brent_search(
                residual_func=residual_kg_s,
                bounds=[low_pressure_pa, high_pressure_pa],
                tolerance=CHAMBER_PRESSURE_RESIDUAL_TOLERANCE_KG_S,
                x_tolerance=CHAMBER_PRESSURE_TOLERANCE_PA,
            )

            converged = (
                abs(new_chamber_pressure_pa - chamber_pressure_pa)
                <= CHAMBER_PRESSURE_TOLERANCE_PA
            )

            chamber_pressure_pa = new_chamber_pressure_pa

            if converged:
                break

        return chamber_pressure_pa




    def calculate_state(
        self,
        ox_mdot_kg_s: float,
        fuel_mdot_kg_s: float,
        ambient_pressure_pa: float,
        previous_state: EngineState,
        chamber_pressure_pa: float | None = None,
   ) -> EngineState:
        """
        Build the engine state for a timestep.

        chamber_pressure_pa should be the value from solve_chamber_pressure_pa,
        computed against the same tank states the mass flows came from. If it is
        omitted the old explicit relation is used, which lags the chamber by one
        timestep and is only stable while the injector pressure drop stays large.
        """

        min_mdot_kg_s = MINIMUM_MDOT_KG_S

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


        mixture_ratio = ox_mdot_kg_s / fuel_mdot_kg_s
        total_mdot_kg_s = ox_mdot_kg_s + fuel_mdot_kg_s

        cea = self.make_cea_obj()

        if chamber_pressure_pa is None:

            # Legacy explicit path. Kept only as a fallback - it lags the chamber
            # by a timestep and goes unstable once the injector pressure drop gets
            # small. Prefer solve_chamber_pressure_pa.
            chamber_pressure_pa = previous_state.chamber_pressure_pa or ATMOSPHERE_PRESSURE_PA

            pc_tolerance = 1e-4

            for _ in range(CHAMBER_PRESSURE_OUTER_ITERATIONS):
                cstar_delivered_m_s = self.cstar_efficiency * cea.get_Cstar(
                    Pc=chamber_pressure_pa,
                    MR=mixture_ratio
                )

                new_chamber_pressure_pa = (
                    total_mdot_kg_s * cstar_delivered_m_s
                    / self.geometry.nozzle_throat_area_m2
                )

                # compare against the previous guess before overwriting it - doing
                # this the other way round made the difference identically zero, so
                # the loop always broke on the first pass
                converged = (
                    abs(new_chamber_pressure_pa - chamber_pressure_pa)
                    / new_chamber_pressure_pa < pc_tolerance
                )

                chamber_pressure_pa = new_chamber_pressure_pa

                if converged:
                    break

        cstar_ideal_m_s = cea.get_Cstar(Pc=chamber_pressure_pa, MR=mixture_ratio)
        cstar_delivered_m_s = self.cstar_efficiency * cstar_ideal_m_s


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
