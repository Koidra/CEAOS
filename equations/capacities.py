import math

from configs import Constants
from data_models import Inputs


def canopy_heat_capacity(inputs: Inputs) -> float:
    """The heat capacity of the canopy
    Equation 8.20

    :param LAI: the leaf area index
    :return: [J K^(-1) m^(-2)]
    """
    return Constants.Global.cap_Leaf * inputs.LAI


def internal_external_canopy_heat_capacity(lumped_cover_heat_capacity: float) -> float:
    """The heat capacity of the canopy
    Equation 8.21

    :param lumped_cover_heat_capacity: the total heat capacity of the lumped cover
    :return: [J K^(-1) m^(-2)]
    """
    return 0.1 * lumped_cover_heat_capacity


def heating_pipe_heat_capacity():
    # Equation 8.22
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    phi_Pipe_i = Constants.Greenhouse.Heating.phi_Pipe_i
    rho_Steel = Constants.Global.rho_Steel
    rho_Water = Constants.Global.rho_Water
    c_pSteel = Constants.Global.c_pSteel
    c_pWater = Constants.Global.c_pWater
    return 0.25 * math.pi * l_Pipe(
        (phi_Pipe_e ** 2 - phi_Pipe_i ** 2) * rho_Steel * c_pSteel + phi_Pipe_i ** 2 * rho_Water * c_pWater)


def remaining_object_heat_capacity(h_obj, rho_obj, c_p_obj):
    # Equation 8.23
    return h_obj * rho_obj * c_p_obj


def air_compartment_water_vapour_capacity(inputs: Inputs):
    # Equation 8.25
    M_Water = Constants.Global.M_Water
    h_Air = Constants.Greenhouse.Construction.h_Air
    R = Constants.Global.R
    air_t = inputs.air_t
    return M_Water * h_Air / ((air_t + 273.15) * R)
