from collections import Callable
from typing import List
from numpy import ndarray as vec
from numpy import ndarray as mat


def rk4(A: mat, x: vec, h: float) -> vec:
    """
    This is an implementation of Runge-Kutta 4th order method
    A: square matrix encapsulating the dynamics dx/dt = A * x. A is assumed to be constant during the step.
    x: state variable at time t
    return x at t+h
    """
    k1 = A * x
    k2 = A * (x + h * k1 / 2)
    k3 = A * (x + h * k2 / 2)
    k4 = A * (x + h * k3)
    return x + h/6 * (k1 + 2*k2 + 2*k3 + k4)


def multistep_solve(f: Callable, xs: List[vec], h: float) -> vec:
    raise NotImplementedError()