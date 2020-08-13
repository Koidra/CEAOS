from constants import *
from crop.tomato.equations.inhibitions import fruit_flow_inhibition
from crop.tomato.equations.utils import fruit_development, smoothed_conditional_function


def fruit_set_of_first_development_stage(last_24_canopy_t, carbohydrate_flow_BufFruits):
    """
    Equation 9.29, B.4
    number_flow_BufFruit_1 = S_number_flow_BufFruit_1_carbohydrate_flow_BufFruits * maximum_fruit_set
    Returns: fruit set of the first development stage [fruits m^-2 s^-1]
    """
    return smoothed_conditional_function(carbohydrate_flow_BufFruits, -58.9, 0.05) * maximum_fruit_set_dependency_on_temperature(last_24_canopy_t)


def maximum_fruit_set_dependency_on_temperature(last_24_canopy_t):
    """
    Equation 9.30
    max_number_flow_BufFruit_1 = PLANT_DENSITY * (MAX_FRUIT_SET_REGRESSION_COEF_1 + MAX_FRUIT_SET_REGRESSION_COEF_2*last_24_canopy_t)
    Returns: maximum fruit set dependency on temperature [fruits m^-2 s^-1]
    """
    return PLANT_DENSITY * (MAX_FRUIT_SET_REGRESSION_COEF_1 + MAX_FRUIT_SET_REGRESSION_COEF_2 * last_24_canopy_t)


def fruit_flow_through_fruit_development_stage(jth: int, number_Fruits, sum_canopy_t, last_24_canopy_t):
    """
    Equation 9.31
    number_flow_Fruit_j_Fruit_jplus = fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * fruit_flow_inhibition_rate * number_fruit_j
    Returns: fruit flow through fruit development stage j [fruits m^-2 s^-1]
    """
    fruit_development_rate = fruit_development(last_24_canopy_t)
    fruit_flow_inhibition_rate = fruit_flow_inhibition(sum_canopy_t)
    number_fruit_j = number_Fruits[jth]
    return fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * fruit_flow_inhibition_rate * number_fruit_j
