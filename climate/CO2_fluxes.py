from ..constants import *
from ..coefficients import Coefficients as coefs
from ..data_models import States, Setpoints, Weather
from .heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air
from .utils import total_side_vents_ventilation_rates, total_roof_ventilation_rates, \
    thermal_screen_air_flux_rate, air_flux


def greenhouse_air_and_above_thermal_screen_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    f_AirTop = thermal_screen_air_flux_rate(setpoints, states, weather)
    return air_flux(f_AirTop, states.air_CO2, states.above_thermal_screen_CO2)


def greenhouse_air_and_outdoor_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = total_side_vent_rate + f_VentForced
    return air_flux(f_AirOut, states.air_CO2, weather.outdoor_CO2)


def above_thermal_screen_and_outdoor_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    f_TopOut = total_roof_ventilation_rates(setpoints, states, weather)
    return air_flux(f_TopOut, states.above_thermal_screen_CO2, weather.outdoor_CO2)


def heat_blower_to_greenhouse_air_CO2_flux(setpoints: Setpoints):
    # Equation 8.54
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return ETA_HEATCO2 * sensible_heat_flux_BlowAir


def external_CO2_added(setpoints: Setpoints):
    # Equation 8.77
    cap_extCO2 = coefs.ActiveClimateControl.cap_extCO2
    floor_area = coefs.Construction.floor_area
    return setpoints.U_ExtCO2 * cap_extCO2 / floor_area
