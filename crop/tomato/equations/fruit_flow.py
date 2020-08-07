from constants import *
from crop.tomato.equations.inhibitions import fruit_flow_inhibition
from crop.tomato.equations.utils import fruit_development


def fruit_set_of_first_development_stage():
    """
    Equation 9.29
    number_flow_BufFruit_1 =
    Returns: fruit set of the first development stage [fruits m^-2 s^-1]
    """
    pass


def maximum_fruit_set_dependency_on_temperature():
    """
    Equation 9.30
    max_number_flow_BufFruit_1 = PLANT_DENSITY * (MAX_FRUIT_SET_REGRESSION_COEF_1 + MAX_FRUIT_SET_REGRESSION_COEF_2*_24_mean_temperature)
    Returns: maximum fruit set dependency on temperature [fruits m^-2 s^-1]
    """
    return PLANT_DENSITY * (MAX_FRUIT_SET_REGRESSION_COEF_1 + MAX_FRUIT_SET_REGRESSION_COEF_2*_24_mean_temperature)


def fruit_flow_through_fruit_development_stage(jth: int):
    """
    Equation 9.31
    number_flow_Fruit_j_Fruit_jplus = fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * fruit_flow_inhibition_rate * number_fruit_j
    Returns: fruit flow through fruit development stage j [fruits m^-2 s^-1]
    """
    fruit_development_rate = fruit_development()
    fruit_flow_inhibition_rate = fruit_flow_inhibition()
    number_fruit_j = 0
    return fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * fruit_flow_inhibition_rate * number_fruit_j