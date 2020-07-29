import math

from configs import Constants
from data_models import Setpoints, Inputs


def air_density():
    # Equation 8.24
    rho_Air0 = Constants.Global.rho_Air0
    g = Constants.Global.g
    M_Air = Constants.Global.M_Air
    h_Elevation = Constants.Greenhouse.Construction.h_Elevation
    R = Constants.Global.R
    return rho_Air0 * math.exp(g * M_Air * h_Elevation / (293.15 * R))


def pad_and_fan_system_ventilation_flux(setpoints: Setpoints):
    # Equation 8.59
    U_Pad = setpoints.U_Pad
    phi_Pad = Constants.Greenhouse.ActiveClimateControl.phi_Pad
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    return U_Pad * phi_Pad / A_Flr


def mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.63
    U_MechCool = setpoints.U_MechCool
    COP_MechCool = Constants.Greenhouse.ActiveClimateControl.COP_MechCool
    P_MechCool = Constants.Greenhouse.ActiveClimateControl.P_MechCool
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    air_t = inputs.air_t
    mechcool_t = inputs.mechcool_t
    delta_H = Constants.Global.delta_H
    VP_Air = saturation_vapour_pressure(air_t)
    VP_MechCool = saturation_vapour_pressure(mechcool_t)
    return (U_MechCool * COP_MechCool * P_MechCool / A_Flr) / (
            air_t - mechcool_t + 6.5E-9 * delta_H * (VP_Air - VP_MechCool))


def roof_ventilation_natural_ventilation_rate(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.65
    U_Roof = setpoints.U_Roof
    max_area_roof_ventilation = Constants.Greenhouse.Ventilation.A_Roof  # TODO: Need to re-check
    C_d = discharge_coefficients(setpoints, 'd')
    C_w = discharge_coefficients(setpoints, 'w')
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    g = Constants.Global.g
    h_Vent = Constants.Greenhouse.Construction.h_Vent
    air_t = inputs.air_t
    outdoor_t = inputs.outdoor_t
    mean_t = (air_t + outdoor_t) / 2
    v_Wind = inputs.v_Wind
    return U_Roof * max_area_roof_ventilation * C_d * math.sqrt(
        g * h_Vent * (air_t - outdoor_t) / (2 * (mean_t + 273.15)) + C_w * v_Wind ** 2) / (2 * A_Flr)


def roof_and_side_vents_ventilation_rate(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.66
    C_d = discharge_coefficients(setpoints, 'd')
    C_w = discharge_coefficients(setpoints, 'w')
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    A_U_Roof = roof_apertures(setpoints)
    A_U_Side = sidewall_vents_apertures(setpoints)
    g = Constants.Global.g
    air_t = inputs.air_t
    outdoor_t = inputs.outdoor_t
    mean_t = (air_t + outdoor_t) / 2
    h_SideRoof = Constants.Greenhouse.Construction.h_SideRoof
    v_Wind = inputs.v_Wind
    return (C_d / A_Flr) * math.sqrt((A_U_Roof * A_U_Side / math.sqrt(A_U_Roof ** 2 + A_U_Side ** 2)) ** 2 * (
            2 * g * h_SideRoof * (air_t - outdoor_t) / (mean_t + 273.15)) + (
                                             (A_U_Roof + A_U_Side) / 2) ** 2 * C_w * v_Wind ** 2)


def sidewall_ventilation_rate(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.67
    C_d = discharge_coefficients(setpoints, 'd')
    C_w = discharge_coefficients(setpoints, 'w')
    A_U_Side = sidewall_vents_apertures(setpoints)
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    v_Wind = inputs.v_Wind
    return C_d * A_U_Side * v_Wind * math.sqrt(C_w) / (2 * A_Flr)


def roof_apertures(setpoints: Setpoints):
    # Equation 8.68
    U_Roof = setpoints.U_Roof
    max_area_roof_ventilation = Constants.Greenhouse.Ventilation.A_Roof  # TODO: Need to re-check
    return U_Roof * max_area_roof_ventilation


def sidewall_vents_apertures(setpoints: Setpoints):
    # Equation 8.69
    U_Side = setpoints.U_Roof
    max_area_sidewall_ventilation = Constants.Greenhouse.Ventilation.A_Side  # TODO: Need to re-check
    return U_Side * max_area_sidewall_ventilation


def ventilation_rate_reduce_factor():
    # Equation 8.70
    sigma_InsScr = Constants.Greenhouse.Ventilation.sigma_InsScr
    return sigma_InsScr * (2 - sigma_InsScr)


def greenhouse_leakage_rate(inputs: Inputs):
    # Equation 8.71
    v_Wind = inputs.v_Wind
    c_leakage = Constants.Greenhouse.Construction.c_leakage
    if v_Wind < 0.25:
        return 0.25 * c_leakage
    else:
        return c_leakage * v_Wind


def total_roof_ventilation_rates(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.72
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    eta_Roof_Thr = Constants.Global.eta_Roof_Thr
    U_ThScr = setpoints.U_ThScr
    eta_InsScr = ventilation_rate_reduce_factor()
    f_mark_VentRoof = roof_ventilation_natural_ventilation_rate(setpoints, inputs)
    f_mark_VentRoofSide = roof_and_side_vents_ventilation_rate(setpoints, inputs)
    f_leakage = greenhouse_leakage_rate(inputs)
    if eta_Roof >= eta_Roof_Thr:
        return eta_InsScr * f_mark_VentRoof + 0.5 * f_leakage
    else:
        return eta_InsScr * (
                U_ThScr * f_mark_VentRoof + (1 - U_ThScr) * f_mark_VentRoofSide * eta_Roof) + 0.5 * f_leakage


def total_side_vents_ventilation_rates(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.73
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    eta_Side = 0  # Note: line 611 / setGlAux / GreenLight
    eta_Roof_Thr = Constants.Global.eta_Roof_Thr
    U_ThScr = setpoints.U_ThScr
    eta_InsScr = ventilation_rate_reduce_factor()
    f_mark_VentSide = sidewall_ventilation_rate(setpoints, inputs)
    f_mark_VentRoofSide = roof_and_side_vents_ventilation_rate(setpoints, inputs)
    f_leakage = greenhouse_leakage_rate(inputs)
    if eta_Roof >= eta_Roof_Thr:
        return eta_InsScr * f_mark_VentSide + 0.5 * f_leakage
    else:
        return eta_InsScr * (U_ThScr * f_mark_VentSide + (1 - U_ThScr) * f_mark_VentRoofSide * eta_Side) + 0.5 * f_leakage


def discharge_coefficients(setpoints: Setpoints, type: str):
    # Equation 8.74, 8.75
    U_ShScr = setpoints.U_ShScr
    if type == 'd':
        C_Gh_d = Constants.Greenhouse.Construction.C_Gh_d
        eta_ShScrC_d = Constants.Greenhouse.Ventilation.eta_ShScrC_d
        return C_Gh_d * (1 - eta_ShScrC_d * U_ShScr)
    else:
        C_Gh_w = Constants.Greenhouse.Construction.C_Gh_w
        eta_ShScrC_w = Constants.Greenhouse.Ventilation.eta_ShScrC_w
        return C_Gh_w * (1 - eta_ShScrC_w * U_ShScr)


def forced_ventilation(setpoints: Setpoints):
    # Equation 8.76
    U_VentForced = setpoints.U_VentForced
    phi_VentForced = Constants.Greenhouse.ActiveClimateControl.phi_VentForced
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    return U_VentForced * phi_VentForced / A_Flr


# def depth_d_soil_temperature(d_Soil, inputs: Inputs):
#     # Equation 8.78
#     soil_mean_t = inputs.soil_mean_t
#     a_0 = inputs.a_0
#     D = damping_depth()
#     omega = Constants.Global.omega
#     t = inputs.t
#     beta = inputs.beta
#     return soil_mean_t + a_0 * math.exp(-d_Soil / D) * math.sin(omega * t - d_Soil / D + beta)


# def damping_depth():
#     # Equation 8.79
#     lambda_soil = Constants.Greenhouse.Soil.lambda_soil
#     rho_Soil =
#     c_pSoil =
#     return math.sqrt(2 * lambda_soil / (rho_Soil * c_pSoil))


# def sky_temperature(inputs: Inputs):
#     # Equation 8.80
#     fr_cloud =
#     epsilon_Sky_Clear = clear_sky_FIR_emission_coefficient(inputs)
#     outdoor_t = inputs.outdoor_t
#     sigma = Constants.Global.sigma
#     return ((1 - fr_cloud) * epsilon_Sky_Clear * (outdoor_t + 273.15) ** 4 + fr_cloud * (
#             (outdoor_t + 372.15) ** 4 - 9 / sigma)) ** 0.25 - 273.15


# def fraction_clouds_during_daytime(inputs: Inputs):
#     # Equation 8.81
#     # TODO: why sum ?
#     I_Glob = inputs.I_Glob
#     I_Sol =
#     return I_Glob / I_Sol


def clear_sky_FIR_emission_coefficient(inputs: Inputs):
    # Equation 8.82
    return 0.53 + 6E-3 * inputs.VP_Out ** 0.5

def saturation_vapour_pressure(temp):
    return 610.78 * math.exp(temp / (temp + 238.3) * 17.2694) # Pascal