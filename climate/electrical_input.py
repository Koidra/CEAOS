from coefficients import Coefficients
from data_models import Setpoints, States


def lamp_electrical_input(setpoints: Setpoints):
    """
    Equation A16 [2]
    Args:
        setpoints:

    Returns: the electrical input to the lamp [W m^-2]
    """
    return Coefficients.Lamp.electrical_capacity_lamp*setpoints.U_Lamp


def inter_lamp_electrical_input(states: States, setpoints: Setpoints):
    raise NotImplemented
