from typing import NamedTuple


class Setpoint(NamedTuple):
    # Lumped Cover Layers setpoint
    # Section 8.4, Equation 8.16, ...
    U_ShadingScreen: float
    U_ShadingScreenPer: float
    U_Roof: float
    U_ThermalScreen: float