import math

from coefficients import Coefficients
from data_models import Setpoints, States, Weather


def air_density():
    # Equation 8.24
    density_Air0 = Coefficients.Outside.density_Air0
    g = Coefficients.Outside.g
    M_Air = Coefficients.Outside.M_Air
    elevation_height = Coefficients.Greenhouse.Construction.elevation_height
    R = Coefficients.Outside.R
    return density_Air0 * math.exp(g * M_Air * elevation_height / (293.15 * R))


def pad_and_fan_system_ventilation_flux(setpoints: Setpoints):
    # Equation 8.59
    U_Pad = setpoints.U_Pad
    phi_Pad = Coefficients.Greenhouse.ActiveClimateControl.phi_Pad
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    return U_Pad * phi_Pad / floor_surface


def mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints: Setpoints, states: States):
    # Equation 8.63
    U_MechCool = setpoints.U_MechCool
    COP_MechCool = Coefficients.Greenhouse.ActiveClimateControl.COP_MechCool
    P_MechCool = Coefficients.Greenhouse.ActiveClimateControl.P_MechCool
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    air_t = states.air_t
    mechcool_t = states.mechcool_t
    evaporation_latent_heat = Coefficients.Outside.evaporation_latent_heat
    VP_Air = saturation_vapor_pressure(air_t)
    VP_MechCool = saturation_vapor_pressure(mechcool_t)
    return (U_MechCool * COP_MechCool * P_MechCool / floor_surface) / (
            air_t - mechcool_t + 6.5E-9 * evaporation_latent_heat * (VP_Air - VP_MechCool))


def roof_ventilation_natural_ventilation_rate(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.65
    U_Roof = setpoints.U_Roof
    max_area_roof_ventilation = Coefficients.Greenhouse.Ventilation.A_Roof  # TODO: Need to re-check
    C_d = discharge_coefficients(setpoints, 'd')
    C_w = discharge_coefficients(setpoints, 'w')
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    g = Coefficients.Outside.g
    h_Vent = Coefficients.Greenhouse.Construction.h_Vent
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    mean_t = (air_t + outdoor_t) / 2
    v_Wind = weather.v_Wind
    return U_Roof * max_area_roof_ventilation * C_d * math.sqrt(
        g * h_Vent * (air_t - outdoor_t) / (2 * (mean_t + 273.15)) + C_w * v_Wind ** 2) / (2 * floor_surface)


def roof_and_side_vents_ventilation_rate(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.66
    C_d = discharge_coefficients(setpoints, 'd')
    C_w = discharge_coefficients(setpoints, 'w')
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    A_U_Roof = roof_apertures(setpoints)
    A_U_Side = sidewall_vents_apertures(setpoints)
    g = Coefficients.Outside.g
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    mean_t = (air_t + outdoor_t) / 2
    h_SideRoof = Coefficients.Greenhouse.Construction.h_SideRoof
    v_Wind = weather.v_Wind
    return (C_d / floor_surface) * math.sqrt((A_U_Roof * A_U_Side / math.sqrt(A_U_Roof ** 2 + A_U_Side ** 2)) ** 2 * (
            2 * g * h_SideRoof * (air_t - outdoor_t) / (mean_t + 273.15)) + (
                                             (A_U_Roof + A_U_Side) / 2) ** 2 * C_w * v_Wind ** 2)


def sidewall_ventilation_rate(setpoints: Setpoints, weather: Weather):
    # Equation 8.67
    C_d = discharge_coefficients(setpoints, 'd')
    C_w = discharge_coefficients(setpoints, 'w')
    A_U_Side = sidewall_vents_apertures(setpoints)
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    v_Wind = weather.v_Wind
    return C_d * A_U_Side * v_Wind * math.sqrt(C_w) / (2 * floor_surface)


def roof_apertures(setpoints: Setpoints):
    # Equation 8.68
    U_Roof = setpoints.U_Roof
    max_area_roof_ventilation = Coefficients.Greenhouse.Ventilation.A_Roof  # TODO: Need to re-check
    return U_Roof * max_area_roof_ventilation


def sidewall_vents_apertures(setpoints: Setpoints):
    # Equation 8.69
    U_Side = setpoints.U_Roof
    max_area_sidewall_ventilation = Coefficients.Greenhouse.Ventilation.A_Side  # TODO: Need to re-check
    return U_Side * max_area_sidewall_ventilation


def ventilation_rate_reduce_factor():
    # Equation 8.70
    sigma_InsScr = Coefficients.Greenhouse.Ventilation.sigma_InsScr
    return sigma_InsScr * (2 - sigma_InsScr)


def greenhouse_leakage_rate(weather: Weather):
    # Equation 8.71
    v_Wind = weather.v_Wind
    leakage = Coefficients.Greenhouse.Construction.leakage
    if v_Wind < 0.25:
        return 0.25 * leakage
    else:
        return leakage * v_Wind


def total_roof_ventilation_rates(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.72
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    eta_Roof_Thr = Coefficients.Outside.eta_Roof_Thr
    U_ThScr = setpoints.U_ThScr
    eta_InsScr = ventilation_rate_reduce_factor()
    f_mark_VentRoof = roof_ventilation_natural_ventilation_rate(setpoints, states, weather)
    f_mark_VentRoofSide = roof_and_side_vents_ventilation_rate(setpoints, states, weather)
    f_leakage = greenhouse_leakage_rate(weather)
    if eta_Roof >= eta_Roof_Thr:
        return eta_InsScr * f_mark_VentRoof + 0.5 * f_leakage
    else:
        return eta_InsScr * (
                U_ThScr * f_mark_VentRoof + (1 - U_ThScr) * f_mark_VentRoofSide * eta_Roof) + 0.5 * f_leakage


def total_side_vents_ventilation_rates(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.73
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    eta_Side = 0  # Note: line 611 / setGlAux / GreenLight
    eta_Roof_Thr = Coefficients.Outside.eta_Roof_Thr
    U_ThScr = setpoints.U_ThScr
    eta_InsScr = ventilation_rate_reduce_factor()
    f_mark_VentSide = sidewall_ventilation_rate(setpoints, weather)
    f_mark_VentRoofSide = roof_and_side_vents_ventilation_rate(setpoints, states, weather)
    f_leakage = greenhouse_leakage_rate(weather)
    if eta_Roof >= eta_Roof_Thr:
        return eta_InsScr * f_mark_VentSide + 0.5 * f_leakage
    else:
        return eta_InsScr * (U_ThScr * f_mark_VentSide + (1 - U_ThScr) * f_mark_VentRoofSide * eta_Side) + 0.5 * f_leakage


def discharge_coefficients(setpoints: Setpoints, type: str):
    # Equation 8.74, 8.75
    U_ShScr = setpoints.U_ShScr
    if type == 'd':
        C_Gh_d = Coefficients.Greenhouse.Construction.C_Gh_d
        eta_ShScrC_d = Coefficients.Greenhouse.Ventilation.eta_ShScrC_d
        return C_Gh_d * (1 - eta_ShScrC_d * U_ShScr)
    else:
        C_Gh_w = Coefficients.Greenhouse.Construction.C_Gh_w
        eta_ShScrC_w = Coefficients.Greenhouse.Ventilation.eta_ShScrC_w
        return C_Gh_w * (1 - eta_ShScrC_w * U_ShScr)


def forced_ventilation(setpoints: Setpoints):
    # Equation 8.76
    U_VentForced = setpoints.U_VentForced
    phi_VentForced = Coefficients.Greenhouse.ActiveClimateControl.phi_VentForced
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    return U_VentForced * phi_VentForced / floor_surface


def clear_sky_FIR_emission_coefficient(weather: Weather):
    # Equation 8.82
    return 0.53 + 6E-3 * weather.VP_Out ** 0.5


def saturation_vapor_pressure(temp):
    return 610.78 * math.exp(temp / (temp + 238.3) * 17.2694)  # Pascal
