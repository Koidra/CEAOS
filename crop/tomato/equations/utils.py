import math

from constants import *


def leaf_area_index(carbohydrate_amount_Leaf):
    """
    Equations 9.5
    LAI = SPECIFIC_LEAF_AREA_INDEX * carbohydrate_amount_Leaf
    Returns: leaf area index [m^2 {leaf} m^-2]
    """
    return SPECIFIC_LEAF_AREA_INDEX * carbohydrate_amount_Leaf


def growth_rate_dependency_to_temperature(_24_canopy_t):
    """
    Equations 9.28
    growth_rate = 0.0047 * _24_canopy_t + 0.06
    Returns: growth rate dependency to temperature [-]
    """
    return 0.0047 * _24_canopy_t + 0.06


def fruit_development(_24_canopy_t):
    """
    Equations 9.32
    fruit_development_rate = FRUIT_DEVELOPMENT_RATE_COEF_1 + FRUIT_DEVELOPMENT_RATE_COEF_2 * _24_canopy_t
    Returns: fruit development rate [s^-1]
    """
    return FRUIT_DEVELOPMENT_RATE_COEF_1 + FRUIT_DEVELOPMENT_RATE_COEF_2 * _24_canopy_t


def fruit_growth(jth: int, _24_canopy_t):
    """
    Equations 9.38
    fruit_growth_rate_j = POTENTIAL_FRUIT_DRY_WEIGHT*math.exp(-math.exp(-curve_steepness*(days_after_fruit_set - fruit_development_time)))
    Returns: fruit growth rate [mg {CH2O} fruit^-1 d^-1]
    """
    fruit_development_rate = fruit_development(_24_canopy_t)
    Fruit_Growth_Period = 1/(fruit_development_rate*86400)
    fruit_development_time = -93.4 + 548.0 * Fruit_Growth_Period
    curve_steepness = 1/(2.44 + 403.0 * fruit_development_time)
    days_after_fruit_set = ((jth-1)+0.5)*Fruit_Growth_Period/FRUIT_DEVELOPMENT_STAGES_NUM
    return POTENTIAL_FRUIT_DRY_WEIGHT*math.exp(-math.exp(-curve_steepness*(days_after_fruit_set - fruit_development_time)))


def sum_carbohydrate_flow_BufFruit_conversion(number_Fruits, _24_canopy_t):
    """
    Equations 9.37
    sum_carbohydrates_flow_BufFruit_conversion_factor = 1/sum(number_fruit_j * fruit_growth_rate_j)
    Returns: fruit growth rate [mg {CH2O} fruit^-1 d^-1]
    """
    total = 1e-10
    for i in range(1, FRUIT_DEVELOPMENT_STAGES_NUM):
        total += number_Fruits[i] * fruit_growth(i, _24_canopy_t)
    return 1/total


def maximum_carbohydrates_stored_in_the_leaves():
    """
    Equations 9.46
    max_carbohydrates_flow_Leaf = LAI_Max / SPECIFIC_LEAF_AREA_INDEX
    Returns: maximum allowed carbohydrates stored in the leaves [mg m^-2]
    """
    return LAI_Max / SPECIFIC_LEAF_AREA_INDEX


def smoothed_conditional_function(state, slope, switch):
    """
    Equation B.1
    Args:
        state:
        slope:
        switch:

    Returns:

    """
    return 1/(1+math.exp(slope*(state-switch)))
