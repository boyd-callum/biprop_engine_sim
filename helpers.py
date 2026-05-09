from typing import Callable


def bisection_search(
        residual_func: Callable, 
        bounds: list, 
        max_iterations: int = 100,
        tolerance: float = 1e-3
        ) -> float:

    low = bounds[0]
    high = bounds[1]

    low_residual = residual_func(low)
    high_residual = residual_func(high)

    if low_residual * high_residual > 0.0:
        raise ValueError(
            f"Bisection bounds do not bracket a root: "
            f"f({low})={low_residual}, f({high})={high_residual}"
        )

    iteration = 0

    while iteration < max_iterations:
        
        mid = 0.5 * (low + high)

        mid_residual = residual_func(mid)
        
        # stop once residual is close enough to zero
        if abs(mid_residual) < tolerance:
            break

        # keep the half-interval that still has the sign change
        if low_residual * mid_residual <= 0.0:
            high = mid
            high_residual = mid_residual
        else:
            low = mid
            low_residual = mid_residual
        
        iteration +=1

        if iteration > 100:
            raise RuntimeError(f"Root search failed to converge after 100 iterations. Final residual of {mid_residual}")

    return 0.5 * (low + high)