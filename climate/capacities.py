import math

from coefficients import Coefficients
from data_models import States
from constants import *


def canopy_heat_capacity(states: States) -> float:
    """The heat capacity of the canopy
    Equation 8.20

    :return: [J K^(-1) m^(-2)]
    """
    return CAP_LEAF * states.leaf_area_index


def internal_external_canopy_heat_capacity(lumped_cover_heat_capacity: float) -> float:
    """The heat capacity of the canopy
    Equation 8.21

    :param lumped_cover_heat_capacity: the total heat capacity of the lumped cover
    :return: [J K^(-1) m^(-2)]
    """
    return 0.1 * lumped_cover_heat_capacity


def heating_pipe_heat_capacity():
    # Equation 8.22
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    phi_internal_pipe = Coefficients.Heating.phi_internal_pipe
    return 0.25 * math.pi * pipe_length * \
           ((phi_external_pipe ** 2 - phi_internal_pipe ** 2) * STEEL_DENSITY * C_PSTEEL + phi_internal_pipe ** 2 * WATER_DENSITY * C_PWATER)


def grow_pipe_heat_capacity():
    return 0.25 * math.pi * Coefficients.GrowPipe.pipe_length \
           * ((Coefficients.GrowPipe.phi_external_pipe**2 - Coefficients.GrowPipe.phi_internal_pipe**2)
              * STEEL_DENSITY * C_PSTEEL + Coefficients.GrowPipe.phi_internal_pipe**2 * WATER_DENSITY * C_PWATER)


def remaining_object_heat_capacity(height_obj, density_obj, heat_cap_obj):
    """
    Equation 8.23
    Args:
        height_obj: mean height of the greenhouse object [m]
        density_obj: density of the greenhouse object [kg m-3]
        heat_cap_obj: specific heat capacity of the object [J K-1 kg-1]

    Returns: [J K^-1 m^-2]
    """
    return height_obj * density_obj * heat_cap_obj


def air_compartment_water_vapor_capacity(states: States):
    # Equation 8.25
    air_height = Coefficients.Construction.air_height
    return M_WATER * air_height / ((states.air_t + 273.15) * M_GAS)
