"""The implementation of Lumped Cover Layers' equations

Based on section 8.4
r_: The reflection coefficient of the layer
t_: The transmission coefficient of the layer

The default model contains four cover layers, i.e.
a movable outdoor shading screen (ShScr),
a semi-permanent shading screen (ShScrPer),
the greenhouse roof (Rf) and
a movable indoor thermal screen (ThScr).
"""
import math

from configs import Constants
from data_models import Setpoints, Inputs


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

def shadingscreen_PAR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_ShScr = setpoints.U_ShScr
    U_ShScrPer = setpoints.U_ShScrPer
    tau_ShScrPAR = Constants.Greenhouse.Shadowscreen.tau_ShScrPAR # line 156 / setGlParams / GreenLight
    rho_ShScrPAR = Constants.Greenhouse.Shadowscreen.rho_ShScrPAR # line 153 / setGlParams / GreenLight
    tau_ShScrPerPAR = Constants.Greenhouse.Whitewash.tau_ShScrPerPAR
    rho_ShScrPerPAR = Constants.Greenhouse.Whitewash.rho_ShScrPerPAR
    return (1 - U_ShScr*(1-tau_ShScrPAR))*(1-U_ShScrPer*(1-tau_ShScrPerPAR)) / (1 - U_ShScr*rho_ShScrPAR*U_ShScrPer*rho_ShScrPerPAR)


def shadingscreen_PAR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_ShScr = setpoints.U_ShScr
    U_ShScrPer = setpoints.U_ShScrPer
    tau_ShScrPAR = Constants.Greenhouse.Shadowscreen.tau_ShScrPAR # line 156 / setGlParams / GreenLight
    rho_ShScrPAR = Constants.Greenhouse.Shadowscreen.rho_ShScrPAR # line 153 / setGlParams / GreenLight
    rho_ShScrPerPAR = Constants.Greenhouse.Whitewash.rho_ShScrPerPAR
    return U_ShScr * rho_ShScrPAR + (1 - U_ShScr * (1-tau_ShScrPAR))**2*U_ShScrPer*rho_ShScrPerPAR/(1-U_ShScr*rho_ShScrPAR*U_ShScrPer*rho_ShScrPerPAR)


def roof_thermal_screen_PAR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    tau_RoofPAR = Constants.Greenhouse.Roof.tau_RfPAR
    rho_RoofPAR = Constants.Greenhouse.Roof.rho_RfPAR
    tau_ThScrPAR = Constants.Greenhouse.Thermalscreen.tau_ThScrPAR
    rho_ThScrPAR = Constants.Greenhouse.Thermalscreen.rho_ThScrPAR
    return (1 - U_Roof*(1-tau_RoofPAR))*(1-U_ThScr*(1-tau_ThScrPAR)) / (1 - U_Roof*rho_RoofPAR*U_ThScr*rho_ThScrPAR)


def roof_thermal_screen_PAR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    tau_RoofPAR = Constants.Greenhouse.Roof.tau_RfPAR
    rho_RoofPAR = Constants.Greenhouse.Roof.rho_RfPAR
    rho_ThScrPAR = Constants.Greenhouse.Thermalscreen.rho_ThScrPAR
    return U_Roof * rho_RoofPAR + (1 - U_Roof * (1-tau_RoofPAR))**2*U_ThScr*rho_ThScrPAR/(1-U_Roof*rho_RoofPAR*U_ThScr*rho_ThScrPAR)


def shadingscreen_NIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_ShScr = setpoints.U_ShScr
    U_ShScrPer = setpoints.U_ShScrPer
    tau_ShScrNIR = Constants.Greenhouse.Shadowscreen.tau_ShScrNIR # line 155 / setGlParams / GreenLight
    rho_ShScrNIR = Constants.Greenhouse.Shadowscreen.rho_ShScrNIR # line 152 / setGlParams / GreenLight
    tau_ShScrPerNIR = Constants.Greenhouse.Whitewash.tau_ShScrPerNIR
    rho_ShScrPerNIR = Constants.Greenhouse.Whitewash.rho_ShScrPerNIR
    return (1 - U_ShScr*(1-tau_ShScrNIR))*(1-U_ShScrPer*(1-tau_ShScrPerNIR)) / (1 - U_ShScr*rho_ShScrNIR*U_ShScrPer*rho_ShScrPerNIR)


def shadingscreen_NIR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_ShScr = setpoints.U_ShScr
    U_ShScrPer = setpoints.U_ShScrPer
    tau_ShScrNIR = Constants.Greenhouse.Shadowscreen.tau_ShScrNIR # line 155 / setGlParams / GreenLight
    rho_ShScrNIR = Constants.Greenhouse.Shadowscreen.rho_ShScrNIR # line 152 / setGlParams / GreenLight
    rho_ShScrPerNIR = Constants.Greenhouse.Whitewash.rho_ShScrPerNIR
    return U_ShScr * rho_ShScrNIR + (1 - U_ShScr * (1-tau_ShScrNIR))**2*U_ShScrPer*rho_ShScrPerNIR/(1-U_ShScr*rho_ShScrNIR*U_ShScrPer*rho_ShScrPerNIR)

def roof_thermal_screen_NIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.16
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    tau_RoofNIR = Constants.Greenhouse.Roof.tau_RfNIR
    rho_RoofNIR = Constants.Greenhouse.Roof.rho_RfNIR
    tau_ThScrNIR = Constants.Greenhouse.Thermalscreen.tau_ThScrNIR
    rho_ThScrNIR = Constants.Greenhouse.Thermalscreen.rho_ThScrNIR
    return (1 - U_Roof*(1-tau_RoofNIR))*(1-U_ThScr*(1-tau_ThScrNIR)) / (1 - U_Roof*rho_RoofNIR*U_ThScr*rho_ThScrNIR)


def roof_thermal_screen_NIR_reflection_coefficient(setpoints: Setpoints):
    # Equation 8.17
    U_Roof = setpoints.U_Roof
    U_ThScr = setpoints.U_ThScr
    tau_RoofNIR = Constants.Greenhouse.Roof.tau_RfNIR
    rho_RoofNIR = Constants.Greenhouse.Roof.rho_RfNIR
    rho_ThScrNIR = Constants.Greenhouse.Thermalscreen.rho_ThScrNIR
    return U_Roof * rho_RoofNIR + (1 - U_Roof * (1-tau_RoofNIR))**2*U_ThScr*rho_ThScrNIR/(1-U_Roof*rho_RoofNIR*U_ThScr*rho_ThScrNIR)


def lumped_cover_heat_capacity(setpoints: Setpoints):
    # Equation 8.18
    psi = Constants.Greenhouse.Construction.psi
    U_ShScrPer = setpoints.U_ShScrPer
    h_ShScrPer = Constants.Greenhouse.Whitewash.h_ShScrPer
    rho_ShScrPer = Constants.Greenhouse.Whitewash.rho_ShScrPer
    c_p_ShScrPer = Constants.Greenhouse.Whitewash.c_p_ShScrPer
    h_Rf = Constants.Greenhouse.Roof.h_Rf
    rho_Rf = Constants.Greenhouse.Roof.rho_Rf
    c_p_Rf = Constants.Greenhouse.Roof.c_p_Rf
    return math.cos(psi)*(U_ShScrPer*h_ShScrPer*rho_ShScrPer*c_p_ShScrPer + h_Rf*rho_Rf*c_p_Rf)


def lumped_cover_conductive_heat_flux(setpoints: Setpoints):
    # Equation 8.19
    h_Rf = Constants.Greenhouse.Roof.h_Rf
    lambda_Rf = Constants.Greenhouse.Roof.lambda_Rf
    U_ShScrPer = setpoints.U_ShScrPer
    h_ShScrPer = Constants.Greenhouse.Whitewash.h_ShScrPer
    lambda_ShScrPer = Constants.Greenhouse.Whitewash.lambda_ShScrPer
    return (h_Rf/lambda_Rf + U_ShScrPer*h_ShScrPer/lambda_ShScrPer)**-1