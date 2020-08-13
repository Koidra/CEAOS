from coefficients import Coefficients
from data_models import Setpoints


def lamp_electrical_input(setpoints: Setpoints):
    """
    Equation A16 [2]
    Args:
        setpoints:

    Returns: the electrical input to the lamp [W m^-2]
    """
    return Coefficients.Lamp.electrical_capacity_lamp*setpoints.U_Lamp


def inter_lamp_electrical_input(setpoints: Setpoints):
    """
    Equation A26 [5]
    Args:
        setpoints:

    Returns: Interlight electrical input [W m^{-2}]
    """
    return Coefficients.Interlight.electrical_capacity_inter_lamp*setpoints.U_IntLamp
