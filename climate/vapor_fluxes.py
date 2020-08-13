from climate.heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air
from climate.utils import *


def differentiable_air_to_obj_vapor_flux(vapor_pressure_1: float, vapor_pressure_2: float, heat_exchange_coef: float):
    """
    Equation 8.44
    Args:
        vapor_pressure_1: the vapor pressure of the air of location 1 [Pa]
        vapor_pressure_2: the saturated vapor pressure of object 2 at its temperature [Pa]
        heat_exchange_coef:  heat exchange coefficient between the air of location 1 to object 2 [W m^-2 K^-1]

    Returns:  the vapor flux from the air to an object by condensation [kg m^-2 s^-1]
    """
    return 6.4E-9 * heat_exchange_coef * (vapor_pressure_1 - vapor_pressure_2) \
           / (1 + math.exp(S_MV12 * (vapor_pressure_1 - vapor_pressure_2)))


def general_vapor_flux(air_flux: float, vapor_pressure_1: float, vapor_pressure_2: float, temp_1: float, temp_2: float):
    """
    Equation 8.45
    Args:
        air_flux: the air flux from location 1 to location 2
        vapor_pressure_1: the vapor pressure of the air of location 1 [Pa]
        vapor_pressure_2: the vapor pressure of the air of location 2 [Pa]
        temp_1: the temperature at location 1 (°C)
        temp_2: the temperature at location 2 (°C)

    Returns: vapour flux from location 1 to location 2 [kg m^-2 s^-1]
    """
    return M_WATER * air_flux * (vapor_pressure_1 / (temp_1 + 273.15) + vapor_pressure_2 / (temp_2 + 273.15)) / M_GAS


def fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints: Setpoints):
    # Equation 8.64
    return setpoints.U_Fog * Coefficients.ActiveClimateControl.cap_Fog / Coefficients.Construction.floor_area


def heat_blower_to_greenhouse_air_vapor_flux(setpoints: Setpoints):
    # Equation 8.55
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return ETA_HEATVAP * sensible_heat_flux_BlowAir


def greenhouse_air_to_thermal_screen_vapor_flux(setpoints: Setpoints, states: ClimateStates):
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(states.t_Air - states.t_ThScr) ** 0.33
    vapor_pressure_ThScr = saturation_vapor_pressure(states.t_ThScr)
    return differentiable_air_to_obj_vapor_flux(states.vapor_pressure_Air, vapor_pressure_ThScr, HEC_AirThScr)


def greenhouse_air_to_above_thermal_screen_vapor_flux(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    f_AirTop = thermal_screen_air_flux_rate(setpoints, states, weather)
    return general_vapor_flux(f_AirTop, states.vapor_pressure_Air, states.vapor_pressure_AboveThScr, states.t_Air,
                              states.t_AboveThScr)


def greenhouse_air_to_outdoor_vapor_flux(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = total_side_vent_rate + f_VentForced
    return general_vapor_flux(f_AirOut, states.vapor_pressure_Air, weather.vapor_pressure_outdoor, states.t_Air,
                              weather.t_Outdoor)


def greenhouse_air_to_mechanical_cooling_vapor_flux(states: ClimateStates, setpoints: Setpoints):
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    vapor_pressure_MechCool = saturation_vapor_pressure(states.t_MechCool)
    return differentiable_air_to_obj_vapor_flux(states.vapor_pressure_Air, vapor_pressure_MechCool, HEC_MechAir)


def above_thermal_screen_to_internal_cover_vapor_flux(states: ClimateStates):
    HEC_TopCov_in = Coefficients.Construction.c_HECin * (states.t_AboveThScr - states.t_Cov_internal) ** 0.33 \
                    * Coefficients.Construction.cover_area / Coefficients.Construction.floor_area
    vapor_pressure_Cov_internal = saturation_vapor_pressure(states.t_Cov_internal)
    return differentiable_air_to_obj_vapor_flux(states.vapor_pressure_AboveThScr, vapor_pressure_Cov_internal,
                                                HEC_TopCov_in)


def above_thermal_screen_to_outdoor_vapor_flux(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    f_TopOut = total_roof_ventilation_rates(setpoints, states, weather)
    return general_vapor_flux(f_TopOut, states.vapor_pressure_AboveThScr, weather.vapor_pressure_outdoor,
                              states.t_AboveThScr, weather.t_Outdoor)
