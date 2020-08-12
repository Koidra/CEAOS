import math

from coefficients import Coefficients
from data_models import Setpoints, States, Weather
from constants import *


def air_density():
    # Equation 8.24
    elevation_height = Coefficients.Construction.elevation_height
    return DENSITY_AIR0 * math.exp(GRAVITY * M_AIR * elevation_height / (293.15 * M_GAS))


def thermal_screen_air_flux_rate(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.41
    U_ThScr = setpoints.U_ThScr
    thScr_flux_coefficient = Coefficients.Thermalscreen.thScr_flux_coefficient
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    density_air = air_density()
    elevation_height = Coefficients.Construction.elevation_height
    pressure = 101325 * (1 - 2.5577e-5 * elevation_height) ** 5.25588
    density_Top = M_AIR * pressure / ((states.above_thermal_screen_t + 273.15) * M_GAS)
    density_Out = density_Top  # = rho_Top, line 715 / setGlAux / GreenLight
    density_mean_Air = (density_air + density_Out) / 2
    return U_ThScr * thScr_flux_coefficient * abs(air_t-outdoor_t) ** 0.66 + \
           (1-U_ThScr) * (0.5 * density_mean_Air * (1 - U_ThScr) * GRAVITY * abs(density_air - density_Out)) ** 0.5 / density_mean_Air


def air_flux(air_flow: float, co2_source: float, co2_target: float) -> float:
    """
    Equation 8.46
    Args:
        co2_source, co2_target: CO2-concentration at location (mg m^-3)
        air_flow: the air flux from location 1 to location 2 (m^3 m^-2 s^-1)
    return: CO2 flux accompanying an air flux  from location 1 to location 2 [mg m^-2 s^-1]
    """
    return air_flow * (co2_source - co2_target)


def mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints: Setpoints, states: States):
    # Equation 8.63
    U_MechCool = setpoints.U_MechCool
    perf_MechCool_coef = Coefficients.ActiveClimateControl.perf_MechCool_coef
    ele_cap_MechCool = Coefficients.ActiveClimateControl.ele_cap_MechCool
    floor_area = Coefficients.Construction.floor_area
    air_t = states.air_t
    mechcool_t = states.mechcool_t
    air_vp = saturation_vapor_pressure(air_t)
    MechCool_vp = saturation_vapor_pressure(mechcool_t)
    return (U_MechCool * perf_MechCool_coef * ele_cap_MechCool / floor_area) / \
           (air_t - mechcool_t + 6.5E-9 * EVAPORATION_LATENT_HEAT * (air_vp - MechCool_vp))


def roof_ventilation_natural_ventilation_rate(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.65
    U_Roof = setpoints.U_Roof
    max_area_roof_ventilation = Coefficients.Ventilation.A_Roof  # TODO: Need to re-check
    discharge_coef = discharge_coefficients(setpoints, 'd')
    global_wind_pressure_coef = discharge_coefficients(setpoints, 'w')
    floor_area = Coefficients.Construction.floor_area
    vent_vertical_dimension = Coefficients.Construction.vent_vertical_dimension
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    mean_t = (air_t + outdoor_t) / 2
    v_Wind = weather.v_Wind
    return U_Roof * max_area_roof_ventilation * discharge_coef * \
           math.sqrt(GRAVITY * vent_vertical_dimension * (air_t - outdoor_t) / (2 * (mean_t + 273.15)) +
                     global_wind_pressure_coef * v_Wind ** 2) / \
           (2 * floor_area)


def roof_and_side_vents_ventilation_rate(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.66
    discharge_coef = discharge_coefficients(setpoints, 'd')
    global_wind_pressure_coef = discharge_coefficients(setpoints, 'w')
    floor_area = Coefficients.Construction.floor_area
    rf_vents = roof_vents_apertures(setpoints)
    side_vents = sidewall_vents_apertures(setpoints)
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    mean_t = (air_t + outdoor_t) / 2
    side_wall_roof_vent_distance = Coefficients.Construction.side_wall_roof_vent_distance
    v_Wind = weather.v_Wind
    return (discharge_coef / floor_area) * \
           math.sqrt((rf_vents * side_vents / math.sqrt(rf_vents ** 2 + side_vents ** 2)) ** 2
                     * (2 * GRAVITY * side_wall_roof_vent_distance * (air_t - outdoor_t) / (mean_t + 273.15)) +
                     ((rf_vents + side_vents) / 2) ** 2 * global_wind_pressure_coef * v_Wind ** 2)


def sidewall_ventilation_rate(setpoints: Setpoints, weather: Weather):
    # Equation 8.67
    discharge_coef = discharge_coefficients(setpoints, 'd')
    global_wind_pressure_coef = discharge_coefficients(setpoints, 'w')
    side_vents = sidewall_vents_apertures(setpoints)
    floor_area = Coefficients.Construction.floor_area
    v_Wind = weather.v_Wind
    return discharge_coef * side_vents * v_Wind * math.sqrt(global_wind_pressure_coef) / (2 * floor_area)


def roof_vents_apertures(setpoints: Setpoints):
    # Equation 8.68
    U_Roof = setpoints.U_Roof
    max_area_roof_ventilation = Coefficients.Ventilation.A_Roof  # TODO: Need to re-check
    return U_Roof * max_area_roof_ventilation


def sidewall_vents_apertures(setpoints: Setpoints):
    # Equation 8.69
    U_Side = setpoints.U_Roof
    max_area_sidewall_ventilation = Coefficients.Ventilation.A_Side  # TODO: Need to re-check
    return U_Side * max_area_sidewall_ventilation


def ventilation_rate_reduce_factor():
    # Equation 8.70
    porosity_InsScr = Coefficients.Ventilation.porosity_InsScr
    return porosity_InsScr * (2 - porosity_InsScr)


def greenhouse_leakage_rate(weather: Weather):
    # Equation 8.71
    v_Wind = weather.v_Wind
    leakage_coef = Coefficients.Construction.leakage_coef
    if v_Wind < 0.25:
        return 0.25 * leakage_coef
    else:
        return leakage_coef * v_Wind


def total_roof_ventilation_rates(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.72
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    U_ThScr = setpoints.U_ThScr
    ventilation_rate_reduced = ventilation_rate_reduce_factor()
    vent_roof_rate = roof_ventilation_natural_ventilation_rate(setpoints, states, weather)
    vent_roof_side_rate = roof_and_side_vents_ventilation_rate(setpoints, states, weather)
    leakage_rate = greenhouse_leakage_rate(weather)
    if eta_Roof >= ETA_ROOF_THR:
        return ventilation_rate_reduced * vent_roof_rate + 0.5 * leakage_rate
    else:
        return ventilation_rate_reduced * (U_ThScr * vent_roof_rate + (1 - U_ThScr) * vent_roof_side_rate * eta_Roof) + 0.5 * leakage_rate


def total_side_vents_ventilation_rates(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.73
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    eta_Side = 0  # Note: line 611 / setGlAux / GreenLight
    U_ThScr = setpoints.U_ThScr
    ventilation_rate_reduced = ventilation_rate_reduce_factor()
    vent_side_rate = sidewall_ventilation_rate(setpoints, weather)
    vent_roof_side_rate = roof_and_side_vents_ventilation_rate(setpoints, states, weather)
    leakage_rate = greenhouse_leakage_rate(weather)
    if eta_Roof >= ETA_ROOF_THR:
        return ventilation_rate_reduced * vent_side_rate + 0.5 * leakage_rate
    else:
        return ventilation_rate_reduced * (U_ThScr * vent_side_rate + (1 - U_ThScr) * vent_roof_side_rate * eta_Side) + 0.5 * leakage_rate


def discharge_coefficients(setpoints: Setpoints, type: str):
    # Equation 8.74, 8.75
    U_ShScr = setpoints.U_ShScr
    if type == 'd':
        C_Gh_d = Coefficients.Construction.C_Gh_d
        eta_ShScrC_d = Coefficients.Ventilation.eta_ShScrC_d
        return C_Gh_d * (1 - eta_ShScrC_d * U_ShScr)
    else:
        C_Gh_w = Coefficients.Construction.C_Gh_w
        eta_ShScrC_w = Coefficients.Ventilation.eta_ShScrC_w
        return C_Gh_w * (1 - eta_ShScrC_w * U_ShScr)


def forced_ventilation(setpoints: Setpoints):
    # Equation 8.76
    U_VentForced = setpoints.U_VentForced
    ventForced_air_flow = Coefficients.ActiveClimateControl.ventForced_air_flow
    floor_area = Coefficients.Construction.floor_area
    return U_VentForced * ventForced_air_flow / floor_area


def clear_sky_FIR_emission_coefficient(weather: Weather):
    # Equation 8.82
    return 0.53 + 6E-3 * weather.outdoor_vp ** 0.5


def saturation_vapor_pressure(temp):
    # Calculation based on
    # http://www.conservationphysics.org/atmcalc/atmoclc2.pdf
    return 610.78 * math.exp(temp / (temp + 238.3) * 17.2694)  # Pascal
