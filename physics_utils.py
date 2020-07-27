from math import exp, log


"""
Units
 - Absolute humidity: g/kg
 - Relative humidity: between 0 and 1
 - Temperature: Celsius
"""

def saturation_ah(t: float) -> float:
    # Ref: https://carnotcycle.wordpress.com/2012/08/04/how-to-convert-relative-humidity-to-absolute-humidity/
    c = 6.112 * 18.02 / 0.08314
    c = c / 1.2  # the absolute humidity in the reference has unit g/m3 (we use g/kg)
    return c * exp(17.67 * t / (t + 243.5)) / (273.15 + t)

# From below, the unit of rh is [0,1], not %
def absolute_humidity(rh: float, t: float) -> float:
    return rh * saturation_ah(t)

def relative_humidity(ah: float, t: float) -> float:
    return ah / saturation_ah(t)

def dewpoint(rh: float, t: float) -> float:
    """
    Calculate the dewpoint from relative_humidity
    Ref: https://www.omnicalculator.com/physics/dew-point#howto
    """
    a, b = 17.62, 243.12
    alpha = log(rh) + a * t / (b + t)
    return b * alpha / (a - alpha)

