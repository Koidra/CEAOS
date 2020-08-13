import math

from coefficients import Coefficients
from data_models import Setpoints, ClimateStates, Weather
from constants import *


def air_density():
    # Equation 8.24
    return DENSITY_AIR0 * math.exp(GRAVITY * M_AIR * Coefficients.Construction.elevation_height / (293.15 * M_GAS))


def thermal_screen_air_flux_rate(setpoints: Setpoints, states: ClimateStates, weather: Weather):
    # Equation 8.41
    density_air = air_density()
    pressure = 101325 * (1 - 2.5577e-5 * Coefficients.Construction.elevation_height) ** 5.25588
    density_Out = M_AIR * pressure / ((states.t_AboveThScr + 273.15) * M_GAS) # = rho_Top, line 715 / setGlAux / GreenLight
    density_mean_Air = (density_air + density_Out) / 2
    return setpoints.U_ThScr * Coefficients.Thermalscreen.thScr_flux_coefficient * abs(states.t_Air - weather.t_Outdoor) ** 0.66 \
           + (1-setpoints.U_ThScr) \
           * (0.5 * density_mean_Air * (1 - setpoints.U_ThScr) * GRAVITY * abs(density_air - density_Out)) ** 0.5 \
           / density_mean_Air


def air_flux(air_flow: float, co2_source: float, co2_target: float) -> float:
    """
    Equation 8.46
    Args:
        co2_source, co2_target: CO2-concentration at location (mg m^-3)
        air_flow: the air flux from location 1 to location 2 (m^3 m^-2 s^-1)
    return: CO2 flux accompanying an air flux  from location 1 to location 2 [mg m^-2 s^-1]
    """
    return air_flow * (co2_source - co2_target)


def mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints: Setpoints, states: ClimateStates):
    # Equation 8.63
    vapor_pressure_MechCool = saturation_vapor_pressure(states.t_MechCool)
    return (setpoints.U_MechCool * Coefficients.ActiveClimateControl.perf_MechCool_coef
            * Coefficients.ActiveClimateControl.ele_cap_MechCool
            / Coefficients.Construction.floor_area) \
           / (states.t_Air - states.t_MechCool + 6.5E-9 * EVAPORATION_LATENT_HEAT * (states.vapor_pressure_Air - vapor_pressure_MechCool))


def roof_ventilation_natural_ventilation_rate(setpoints: Setpoints, states: ClimateStates, weather: Weather):
    # Equation 8.65
    max_area_roof_ventilation = Coefficients.Ventilation.A_Roof  # TODO: Need to re-check
    discharge_coef = discharge_coefficients(setpoints, 'd')
    global_wind_pressure_coef = discharge_coefficients(setpoints, 'w')
    vent_vertical_dimension = Coefficients.Construction.vent_vertical_dimension
    mean_t = (states.t_Air + weather.t_Outdoor) / 2
    return setpoints.U_Roof * max_area_roof_ventilation * discharge_coef \
           * math.sqrt(GRAVITY * vent_vertical_dimension * (states.t_Air - weather.t_Outdoor) / (2 * (mean_t + 273.15))
                       + global_wind_pressure_coef * weather.v_Wind ** 2) \
           / (2 * Coefficients.Construction.floor_area)


def roof_and_side_vents_ventilation_rate(setpoints: Setpoints, states: ClimateStates, weather: Weather):
    # Equation 8.66
    discharge_coef = discharge_coefficients(setpoints, 'd')
    global_wind_pressure_coef = discharge_coefficients(setpoints, 'w')
    rf_vents = roof_vents_apertures(setpoints)
    side_vents = sidewall_vents_apertures(setpoints)
    mean_t = (states.t_Air + weather.t_Outdoor) / 2
    return (discharge_coef / Coefficients.Construction.floor_area) * \
           math.sqrt((rf_vents * side_vents / math.sqrt(rf_vents ** 2 + side_vents ** 2)) ** 2
                     * (2 * GRAVITY * Coefficients.Construction.side_wall_roof_vent_distance
                        * (states.t_Air - weather.t_Outdoor) / (mean_t + 273.15))
                     + ((rf_vents + side_vents) / 2) ** 2 * global_wind_pressure_coef * weather.v_Wind ** 2)


def sidewall_ventilation_rate(setpoints: Setpoints, weather: Weather):
    # Equation 8.67
    discharge_coef = discharge_coefficients(setpoints, 'd')
    global_wind_pressure_coef = discharge_coefficients(setpoints, 'w')
    side_vents = sidewall_vents_apertures(setpoints)
    return discharge_coef * side_vents * weather.v_Wind * math.sqrt(global_wind_pressure_coef) \
           / (2 * Coefficients.Construction.floor_area)


def roof_vents_apertures(setpoints: Setpoints):
    # Equation 8.68
    max_area_roof_ventilation = Coefficients.Ventilation.A_Roof  # TODO: Need to re-check
    return setpoints.U_Roof * max_area_roof_ventilation


def sidewall_vents_apertures(setpoints: Setpoints):
    # Equation 8.69
    U_Side = setpoints.U_Roof
    max_area_sidewall_ventilation = Coefficients.Ventilation.A_Side  # TODO: Need to re-check
    return U_Side * max_area_sidewall_ventilation


def ventilation_rate_reduce_factor():
    # Equation 8.70
    return Coefficients.Ventilation.porosity_InsScr * (2 - Coefficients.Ventilation.porosity_InsScr)


def greenhouse_leakage_rate(weather: Weather):
    # Equation 8.71
    if weather.v_Wind < 0.25:
        return 0.25 * Coefficients.Construction.leakage_coef
    else:
        return Coefficients.Construction.leakage_coef * weather.v_Wind


def total_roof_ventilation_rates(setpoints: Setpoints, states: ClimateStates, weather: Weather):
    # Equation 8.72
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    ventilation_rate_reduced = ventilation_rate_reduce_factor()
    vent_roof_rate = roof_ventilation_natural_ventilation_rate(setpoints, states, weather)
    vent_roof_side_rate = roof_and_side_vents_ventilation_rate(setpoints, states, weather)
    leakage_rate = greenhouse_leakage_rate(weather)
    if eta_Roof >= ETA_ROOF_THR:
        return ventilation_rate_reduced * vent_roof_rate + 0.5 * leakage_rate
    else:
        return ventilation_rate_reduced * (setpoints.U_ThScr * vent_roof_rate
                                           + (1 - setpoints.U_ThScr) * vent_roof_side_rate * eta_Roof) \
               + 0.5 * leakage_rate


def total_side_vents_ventilation_rates(setpoints: Setpoints, states: ClimateStates, weather: Weather):
    # Equation 8.73
    eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
    eta_Side = 0  # Note: line 611 / setGlAux / GreenLight
    ventilation_rate_reduced = ventilation_rate_reduce_factor()
    vent_side_rate = sidewall_ventilation_rate(setpoints, weather)
    vent_roof_side_rate = roof_and_side_vents_ventilation_rate(setpoints, states, weather)
    leakage_rate = greenhouse_leakage_rate(weather)
    if eta_Roof >= ETA_ROOF_THR:
        return ventilation_rate_reduced * vent_side_rate + 0.5 * leakage_rate
    else:
        return ventilation_rate_reduced * (setpoints.U_ThScr * vent_side_rate
                                           + (1 - setpoints.U_ThScr) * vent_roof_side_rate * eta_Side) \
               + 0.5 * leakage_rate


def discharge_coefficients(setpoints: Setpoints, type: str):
    # Equation 8.74, 8.75
    if type == 'd':
        return Coefficients.Construction.C_Gh_d * (1 - Coefficients.Ventilation.eta_ShScrC_d * setpoints.U_ShScr)
    else:
        return Coefficients.Construction.C_Gh_w * (1 - Coefficients.Ventilation.eta_ShScrC_w * setpoints.U_ShScr)


def forced_ventilation(setpoints: Setpoints):
    # Equation 8.76
    return setpoints.U_VentForced * Coefficients.ActiveClimateControl.ventForced_air_flow / Coefficients.Construction.floor_area


def clear_sky_FIR_emission_coefficient(weather: Weather):
    # Equation 8.82
    return 0.53 + 6E-3 * weather.vapor_pressure_outdoor ** 0.5


def saturation_vapor_pressure(temp):
    # Calculation based on
    # http://www.conservationphysics.org/atmcalc/atmoclc2.pdf
    return 610.78 * math.exp(temp / (temp + 238.3) * 17.2694)  # Pascal
