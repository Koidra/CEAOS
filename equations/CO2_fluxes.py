from coefficients import Coefficients
from data_models import States, Setpoints, Weather
from equations.utils import total_side_vents_ventilation_rates, total_roof_ventilation_rates
from equations.heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air, \
    thermal_screen_air_flux_rate


def general_CO2_flux(f_12: float, CO2_1, CO2_2):
    # Equation 8.46
    return f_12 * (CO2_1 - CO2_2)


# Not included in GreenLight
# def pad_fan_system_and_greenhouse_air_CO2_flux(states: States):
#     f_PadAir =
#     CO2_Pad =
#     CO2_Air = states.CO2_Air
#     return general_CO2_flux(f_PadAir, CO2_Pad, CO2_Air)


def greenhouse_air_and_above_thermal_screen_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    f_ThScr = thermal_screen_air_flux_rate(setpoints, states, weather)
    f_AirTop = f_ThScr
    CO2_Air = states.CO2_Air
    CO2_Top = states.CO2_Top
    return general_CO2_flux(f_AirTop, CO2_Air, CO2_Top)


def greenhouse_air_and_outdoor_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    f_VentSide = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = f_VentSide + f_VentForced
    CO2_Air = states.CO2_Air
    CO2_Out = weather.CO2_Out
    return general_CO2_flux(f_AirOut, CO2_Air, CO2_Out)


def above_thermal_screen_and_outdoor_CO2_flux(states: States, setpoints: Setpoints, weather: Weather):
    f_VentRoof = total_roof_ventilation_rates(setpoints, states, weather)
    f_TopOut = f_VentRoof
    CO2_Top = states.CO2_Top
    CO2_Out = weather.CO2_Out
    return general_CO2_flux(f_TopOut, CO2_Top, CO2_Out)


def heat_blower_to_greenhouse_air_CO2_flux(setpoints: Setpoints):
    # Equation 8.54
    eta_HeatCO2 = Coefficients.Outside.eta_HeatCO2
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return eta_HeatCO2 * sensible_heat_flux_BlowAir


def external_CO2_added(setpoints: Setpoints):
    # Equation 8.77
    U_ExtCO2 = setpoints.U_ExtCO2
    extCO2_capacity = Coefficients.Greenhouse.ActiveClimateControl.extCO2_capacity
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    return U_ExtCO2 * extCO2_capacity / floor_surface
