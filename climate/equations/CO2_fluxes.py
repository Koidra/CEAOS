from coefficients import Coefficients
from constants import *
from data_models import States, Setpoints, Weather
from climate.equations.utils import total_side_vents_ventilation_rates, total_roof_ventilation_rates, \
    thermal_screen_air_flux_rate
from climate.equations.heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air


def general_CO2_flux(f_12: float, CO2_location_1: float, CO2_location_2: float) -> float:
    """
    Equation 8.46
    return:  CO2 flux accompanying an air flux [mg m^-2 s^-1]
    """
    return f_12 * (CO2_location_1 - CO2_location_2)


def greenhouse_air_and_above_thermal_screen_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    thScr_air_flux_rate = thermal_screen_air_flux_rate(setpoints, states, weather)
    f_AirTop = thScr_air_flux_rate
    CO2_air = states.CO2_air
    CO2_above_thermal_screen = states.CO2_above_thermal_screen
    return general_CO2_flux(f_AirTop, CO2_air, CO2_above_thermal_screen)


def greenhouse_air_and_outdoor_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = total_side_vent_rate + f_VentForced
    CO2_air = states.CO2_air
    CO2_Out = weather.CO2_Out
    return general_CO2_flux(f_AirOut, CO2_air, CO2_Out)


def above_thermal_screen_and_outdoor_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    total_roof_vent_rate = total_roof_ventilation_rates(setpoints, states, weather)
    f_TopOut = total_roof_vent_rate
    CO2_above_thermal_screen = states.CO2_above_thermal_screen
    CO2_Out = weather.CO2_Out
    return general_CO2_flux(f_TopOut, CO2_above_thermal_screen, CO2_Out)


def heat_blower_to_greenhouse_air_CO2_flux(setpoints: Setpoints):
    # Equation 8.54
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return ETA_HEATCO2 * sensible_heat_flux_BlowAir


def external_CO2_added(setpoints: Setpoints):
    # Equation 8.77
    U_ExtCO2 = setpoints.U_ExtCO2
    cap_extCO2 = Coefficients.ActiveClimateControl.cap_extCO2
    floor_area = Coefficients.Construction.floor_area
    return U_ExtCO2 * cap_extCO2 / floor_area
