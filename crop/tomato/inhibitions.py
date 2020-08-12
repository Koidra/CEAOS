import math

from constants import *
from crop.tomato.equations.utils import smoothed_conditional_function


def carbohydrates_saturation_photosynthesis_rate_inhibition(carbohydrate_amount_Buf):
    """
    Equation 9.11
    Returns:  photosynthesis inhibition by carbohydrates saturation rate[-]
    """
    return smoothed_conditional_function(carbohydrate_amount_Buf, 5e-4, 20e3)


def non_optimal_instantaneous_temperature_inhibition(canopy_t):
    """
    Equation B.2
    Returns: growth inhibition by non-optimal instantaneous temperature [-]
    """
    return smoothed_conditional_function(canopy_t, -0.8690, 10) * smoothed_conditional_function(canopy_t, 0.5793, 34)


def non_optimal_24_hour_canopy_temperatures_inhibition(last_24_canopy_t):
    """
    Equation B.3
    Returns: growth inhibition by non-optimal 24-hour mean temperature [-]
    """
    return smoothed_conditional_function(last_24_canopy_t, -1.1587, 15) \
           * smoothed_conditional_function(last_24_canopy_t, 1.3904, 24.5)


def crop_development_stage_inhibition(sum_canopy_t):
    """
    Equation 9.27, B.6
    Returns: The gradual increase in fruit growth rate depending on tomato development stage [-]
    """
    return 0.5 * ((sum_canopy_t / SUM_END_T) + math.sqrt((sum_canopy_t / SUM_END_T) ** 2 + 1e-4)) \
           - 0.5 * ((sum_canopy_t - SUM_END_T / SUM_END_T)
                    + math.sqrt((sum_canopy_t - SUM_END_T / SUM_END_T) ** 2 + 1e-4))


def fruit_flow_inhibition(sum_canopy_t):
    """
    Equation 9.33
    To assure that fruits stay in the first development stage at vegetative stage
    Returns: fruit flow inhibition [-]
    """
    return smoothed_conditional_function(sum_canopy_t, -5e-2, 0)