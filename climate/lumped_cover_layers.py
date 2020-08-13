"""The implementation of Lumped Cover Layers' equations

Based on section 8.4
r_: The reflection coefficient of the layer
t_: The transmission coefficient of the layer

The default model contains three cover layers, i.e.
    - A movable outdoor shading screen (ShScr)
    - The greenhouse roof (Rf)
    - A movable indoor thermal screen (ThScr).
"""
import math

from coefficients import Coefficients
from data_models import Setpoints


def double_layer_cover_transmission_coefficient(transmission_coef_1, transmission_coef_2, reflection_coef_1, reflection_coef_2) -> float:
    """
    The transmission coefficient of a double layer cover
    Equation 8.14
    :param float transmission_coef_1: the transmission coefficients of the first layer
    :param float transmission_coef_2: the transmission coefficients of the second layer
    :param float reflection_coef_1: the reflection coefficients of the first layer
    :param float reflection_coef_2: the reflection coefficients of the second layer
    :return: The transmission coefficient
    """
    return (transmission_coef_1 * transmission_coef_2) / (1 - reflection_coef_1 * reflection_coef_2)


def double_layer_cover_reflection_coefficient(transmission_coef_1, reflection_coef_1, reflection_coef_2) -> float:
    """
    The reflection coefficient
    Equation 8.15
    :param float transmission_coef_1: the transmission coefficients of the first layer
    :param float reflection_coef_1: the reflection coefficients of the first layer
    :param float reflection_coef_2: the reflection coefficients of the second layer
    :return: The reflection coefficient
    """

    return reflection_coef_1 + (transmission_coef_1 * transmission_coef_1 * reflection_coef_2) / (1 - reflection_coef_1 * reflection_coef_2)


def absorption_coefficient(transmission_coef, reflection_coef):
    """The absorption coefficient

    :param transmission_coef: the transmission coefficient
    :param reflection_coef: the reflection coefficient
    :return: The absorption coefficient
    """
    return 1 - (transmission_coef + reflection_coef)

# TODO: refactor this


def roof_thermal_screen_PAR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_PAR_transmission_coef = Coefficients.Roof.roof_PAR_transmission_coefficient
    roof_PAR_reflection_coef = Coefficients.Roof.roof_PAR_reflection_coefficient
    thScr_PAR_transmission_coef = Coefficients.Thermalscreen.thScr_PAR_transmission_coefficient
    thScr_PAR_reflection_coef = Coefficients.Thermalscreen.thScr_PAR_reflection_coefficient
    return (1 - U_Roof * (1 - roof_PAR_transmission_coef)) * \
           (1 - U_ThScr * (1 - thScr_PAR_transmission_coef)) / \
           (1 - U_Roof * roof_PAR_reflection_coef * U_ThScr * thScr_PAR_reflection_coef)


def roof_thermal_screen_PAR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_PAR_transmission_coef = Coefficients.Roof.roof_PAR_transmission_coefficient
    roof_PAR_reflection_coef = Coefficients.Roof.roof_PAR_reflection_coefficient
    thScr_PAR_reflection_coef = Coefficients.Thermalscreen.thScr_PAR_reflection_coefficient
    return U_Roof * roof_PAR_reflection_coef + \
           (1 - U_Roof * (1 - roof_PAR_transmission_coef)) ** 2 * U_ThScr * thScr_PAR_reflection_coef / \
           (1 - U_Roof * roof_PAR_reflection_coef * U_ThScr * thScr_PAR_reflection_coef)


def roof_thermal_screen_NIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_NIR_transmission_coef = Coefficients.Roof.roof_NIR_transmission_coefficient
    roof_NIR_reflection_coef = Coefficients.Roof.roof_NIR_reflection_coefficient
    thScr_NIR_transmission_coef = Coefficients.Thermalscreen.thScr_NIR_transmission_coefficient
    thScr_NIR_reflection_coef = Coefficients.Thermalscreen.thScr_NIR_reflection_coefficient
    return (1 - U_Roof * (1 - roof_NIR_transmission_coef)) *\
           (1 - U_ThScr * (1 - thScr_NIR_transmission_coef)) / \
           (1 - U_Roof * roof_NIR_reflection_coef * U_ThScr * thScr_NIR_reflection_coef)


def roof_thermal_screen_NIR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_NIR_transmission_coef = Coefficients.Roof.roof_NIR_transmission_coefficient
    roof_NIR_reflection_coef = Coefficients.Roof.roof_NIR_reflection_coefficient
    thScr_NIR_reflection_coef = Coefficients.Thermalscreen.thScr_NIR_reflection_coefficient
    return U_Roof * roof_NIR_reflection_coef + \
           (1 - U_Roof * (1 - roof_NIR_transmission_coef)) ** 2 * U_ThScr * thScr_NIR_reflection_coef / \
           (1 - U_Roof * roof_NIR_reflection_coef * U_ThScr * thScr_NIR_reflection_coef)


def lumped_cover_heat_capacity(setpoints: Setpoints):
    # Equation 8.18
    mean_greenhouse_cover_slope = Coefficients.Construction.mean_greenhouse_cover_slope
    roof_thickness = Coefficients.Roof.roof_thickness
    roof_density = Coefficients.Roof.roof_density
    c_p_Rf = Coefficients.Roof.c_p_Rf
    return math.cos(mean_greenhouse_cover_slope)*(roof_thickness*roof_density*c_p_Rf)


def lumped_cover_conductive_heat_flux(setpoints: Setpoints):
    # Equation 8.19
    roof_thickness = Coefficients.Roof.roof_thickness
    roof_heat_conductivity = Coefficients.Roof.roof_heat_conductivity
    return (roof_thickness/roof_heat_conductivity)**-1
