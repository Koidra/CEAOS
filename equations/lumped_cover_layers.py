"""The implementation of Lumped Cover Layers' equations

Based on section 8.4
r_: The reflection coefficient of the layer
t_: The transmission coefficient of the layer

The default model contains three cover layers, i.e.
a movable outdoor shading screen (ShScr),
a semi-permanent shading screen (ShScrPer),
the greenhouse roof (Rf) and
a movable indoor thermal screen (ThScr).
"""
import math

from coefficients import Coefficients
from data_models import Setpoints


def double_layer_cover_transmission_coefficient(tau_1, tau_2, rho_1, rho_2) -> float:
    """
    The transmission coefficient of a double layer cover
    Equation 8.14
    :param float tau_1: the transmission coefficients of the first layer
    :param float tau_2: the transmission coefficients of the second layer
    :param float rho_1: the reflection coefficients of the first layer
    :param float rho_2: the reflection coefficients of the second layer
    :return: The transmission coefficient
    """
    return (tau_1 * tau_2) / (1 - rho_1 * rho_2)


def double_layer_cover_reflection_coefficient(tau_1, rho_1, rho_2) -> float:
    """
    The reflection coefficient
    Equation 8.15
    :param float tau_1: the transmission coefficients of the first layer
    :param float rho_1: the reflection coefficients of the first layer
    :param float rho_2: the reflection coefficients of the second layer
    :return: The reflection coefficient
    """

    return rho_1 + (tau_1 * tau_1 * rho_2) / (1 - rho_1 * rho_2)


def absorption_coefficient(t, r):
    """The absorption coefficient

    :param t: the transmission coefficient
    :param r: the reflection coefficient
    :return: The absorption coefficient
    """
    return 1 - (t + r)

# TODO: refactor this


def roof_thermal_screen_PAR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_PAR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_PAR_transmission_coefficient
    roof_PAR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_PAR_reflection_coefficient
    thScr_PAR_transmission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_PAR_transmission_coefficient
    thScr_PAR_reflection_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_PAR_reflection_coefficient
    return (1 - U_Roof*(1-roof_PAR_transmission_coefficient))*(1-U_ThScr*(1-thScr_PAR_transmission_coefficient)) / (1 - U_Roof*roof_PAR_reflection_coefficient*U_ThScr*thScr_PAR_reflection_coefficient)


def roof_thermal_screen_PAR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_PAR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_PAR_transmission_coefficient
    roof_PAR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_PAR_reflection_coefficient
    thScr_PAR_reflection_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_PAR_reflection_coefficient
    return U_Roof * roof_PAR_reflection_coefficient + (1 - U_Roof * (1-roof_PAR_transmission_coefficient))**2*U_ThScr*thScr_PAR_reflection_coefficient/(1-U_Roof*roof_PAR_reflection_coefficient*U_ThScr*thScr_PAR_reflection_coefficient)


def shadingscreen_NIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_ShScr = setpoints.U_ShScr
    shScr_NIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight
    return (1 - U_ShScr*(1-shScr_NIR_transmission_coefficient)) / (1 - U_ShScr*shScr_NIR_reflection_coefficient)


def shadingscreen_NIR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_ShScr = setpoints.U_ShScr
    shScr_NIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight
    return U_ShScr * shScr_NIR_reflection_coefficient + (1 - U_ShScr * (1-shScr_NIR_transmission_coefficient))**2/(1-U_ShScr*shScr_NIR_reflection_coefficient)


def roof_thermal_screen_NIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_NIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_NIR_transmission_coefficient
    roof_NIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_NIR_reflection_coefficient
    thScr_NIR_transmission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_NIR_transmission_coefficient
    thScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_NIR_reflection_coefficient
    return (1 - U_Roof*(1-roof_NIR_transmission_coefficient))*(1-U_ThScr*(1-thScr_NIR_transmission_coefficient)) / (1 - U_Roof*roof_NIR_reflection_coefficient*U_ThScr*thScr_NIR_reflection_coefficient)


def roof_thermal_screen_NIR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    roof_NIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_NIR_transmission_coefficient
    roof_NIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_NIR_reflection_coefficient
    thScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_NIR_reflection_coefficient
    return U_Roof * roof_NIR_reflection_coefficient + (1 - U_Roof * (1-roof_NIR_transmission_coefficient))**2*U_ThScr*thScr_NIR_reflection_coefficient/(1-U_Roof*roof_NIR_reflection_coefficient*U_ThScr*thScr_NIR_reflection_coefficient)


def lumped_cover_heat_capacity(setpoints: Setpoints):
    # Equation 8.18
    mean_greenhouse_cover_slope = Coefficients.Greenhouse.Construction.mean_greenhouse_cover_slope
    roof_thickness = Coefficients.Greenhouse.Roof.roof_thickness
    roof_density = Coefficients.Greenhouse.Roof.roof_density
    c_p_Rf = Coefficients.Greenhouse.Roof.c_p_Rf
    return math.cos(mean_greenhouse_cover_slope)*(roof_thickness*roof_density*c_p_Rf)


def lumped_cover_conductive_heat_flux(setpoints: Setpoints):
    # Equation 8.19
    roof_thickness = Coefficients.Greenhouse.Roof.roof_thickness
    roof_heat_conductivity = Coefficients.Greenhouse.Roof.roof_heat_conductivity
    return (roof_thickness/roof_heat_conductivity)**-1
