from collections import Callable
from typing import List
import numpy as np
from numpy import ndarray as vec, mat

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


def adams_moulton(As: List[mat], xs: List[vec], h: float) -> vec:
    """
    Adams-Moulton is a linear multistep method which uses the information from the previous steps
    https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Moulton_methods
    """
    a = np.eye(len(xs[0])) - h * 251/720 * As[4]
    b = xs[3] + h * (646 * As[3] * xs[3]
                 - 264 * As[2] * xs[2]
                 + 106 * As[1] * xs[1]
                 - 19 * As[0] * xs[0]) / 720
    return np.linalg.solve(a, b)