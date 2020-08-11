import math

from coefficients import Coefficients
from constants import *
from crop.tomato.equations.utils import leaf_area_index


def PAR_absorbed_by_canopy(carbohydrate_amount_Leaf, outdoor_global_rad):
    """
    Equation 9.17
    canopy_PAR_absorbed = greenhouse_cover_PAR_transmitted + greenhouse_floor_PAR_reflected
    Returns: PAR absorbed by canopy [µmol {photons} m^-2 s^-1]
    """
    greenhouse_cover_PAR_transmitted = PAR_transmitted_by_greenhouse_cover(carbohydrate_amount_Leaf, outdoor_global_rad)
    greenhouse_floor_PAR_reflected = PAR_reflected_by_greenhouse_floor(carbohydrate_amount_Leaf, outdoor_global_rad)
    return greenhouse_cover_PAR_transmitted + greenhouse_floor_PAR_reflected


def PAR_transmitted_by_greenhouse_cover(carbohydrate_amount_Leaf, outdoor_global_rad):
    """
    Equation 9.18
    greenhouse_cover_PAR_transmitted = above_canopy_PAR * (1 - CANOPY_PAR_REFLECTION_COEF) * (1 - math.exp(-CANOPY_PAR_EXTINCTION_COEF*leaf_area_index))
    Returns: PAR transmitted by greenhouse cover [µmol {photons} m^-2 s^-1]
    """
    above_canopy_PAR = PAR_above_the_canopy(outdoor_global_rad)
    return above_canopy_PAR * (1 - CANOPY_PAR_REFLECTION_COEF) * (1 - math.exp(-CANOPY_PAR_EXTINCTION_COEF*leaf_area_index(carbohydrate_amount_Leaf)))


def PAR_above_the_canopy(outdoor_global_rad):
    """
    Equation 9.19
    above_canopy_PAR = COVER_LIGHT_TRANSMISSION * GLOBAL_RADIATION_TO_PAR_CONVERSION * outdoor_global_rad
    Returns: PAR above the canopy [µmol {photons} m^-2 s^-1]
    """
    return COVER_LIGHT_TRANSMISSION * GLOBAL_RADIATION_TO_PAR_CONVERSION * outdoor_global_rad


def PAR_reflected_by_greenhouse_floor(carbohydrate_amount_Leaf, outdoor_global_rad):
    """
    Equation 9.20
    greenhouse_floor_PAR_reflected = floor_PAR_reflection_coef * above_canopy_PAR * (1 - CANOPY_PAR_REFLECTION_COEF) \
                                    * math.exp(-CANOPY_PAR_EXTINCTION_COEF*leaf_area_index) * (1 - math.exp(-FLOOR_PAR_EXTINCTION_COEF*leaf_area_index))
    Returns: PAR reflected by greenhouse floor [µmol {photons} m^-2 s^-1]
    """
    floor_PAR_reflection_coef = Coefficients.Floor.floor_PAR_reflection_coefficient
    above_canopy_PAR = PAR_above_the_canopy(outdoor_global_rad)
    return floor_PAR_reflection_coef * above_canopy_PAR * (1 - CANOPY_PAR_REFLECTION_COEF) \
        * math.exp(-CANOPY_PAR_EXTINCTION_COEF*leaf_area_index(carbohydrate_amount_Leaf)) * (1 - math.exp(-FLOOR_PAR_EXTINCTION_COEF*leaf_area_index(carbohydrate_amount_Leaf)))
