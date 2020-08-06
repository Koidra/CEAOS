from math import exp, log


"""
Units
 - Absolute humidity: g/m3
 - Relative humidity: between 0 and 1
 - Temperature: Celsius
"""

GAS_CONSTANT = 8.314463  # J/K/mol
WATER_MOLECULAR_WEIGHT = 18.02  # g/mol

def _kelvin(t_celsius: float) -> float:
    return t_celsius + 273.15

def saturation_vapor_pressure(t: float) -> float:
    return 0.6112 * exp(17.67 * t / (t + 243.5))

def vapor_pressure(rh: float, t: float) -> float:
    return rh * saturation_vapor_pressure(t)

# From below, the unit of rh is [0,1], not %
def absolute_humidity(rh: float, t: float) -> float:
    n = vapor_pressure(rh, t) / (GAS_CONSTANT * _kelvin(t))  # ideal gas equation: PV = nRT; V=1
    return n * WATER_MOLECULAR_WEIGHT

def relative_humidity(vp: float, t: float) -> float:
    return vp / saturation_vapor_pressure(t)

def dewpoint(rh: float, t: float) -> float:
    """
    Calculate the dewpoint from relative_humidity
    Ref: https://www.omnicalculator.com/physics/dew-point#howto
    """
    a, b = 17.62, 243.12
    alpha = log(rh) + a * t / (b + t)
    return b * alpha / (a - alpha)

