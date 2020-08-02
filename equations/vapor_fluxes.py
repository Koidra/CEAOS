from equations.heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air
from equations.utils import *


def air_to_obj_vapor_flux(VP_1: float, VP_2: float, HEC: float):
    # Equation 8.43
    return 6.4E-9 * HEC * max(VP_1 - VP_2, 0)


def differentiable_air_to_obj_vapor_flux(s_MV12: float, VP_1: float, VP_2: float, HEC: float):
    # Equation 8.44
    return 6.4E-9 * HEC * (VP_1 - VP_2) / (1 + math.exp(s_MV12 * (VP_1 - VP_2)))


def general_vapor_flux(f_12: float, VP_1: float, VP_2: float, object_1_t: float, object_2_t: float):
    # Equation 8.45
    M_Water = Constants.M_Water
    M_Gas = Constants.M_Gas
    return M_Water * f_12 * (VP_1 / (object_1_t + 273.15) + VP_2 / (object_2_t + 273.15)) / M_Gas


def fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints: Setpoints):
    # Equation 8.64
    U_Fog = setpoints.U_Fog
    cap_Fog = Coefficients.ActiveClimateControl.cap_Fog
    floor_area = Coefficients.Construction.floor_area
    return U_Fog * cap_Fog / floor_area


def heat_blower_to_greenhouse_air_vapor_flux(setpoints: Setpoints):
    # Equation 8.55
    eta_HeatVap = Constants.eta_HeatVap
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return eta_HeatVap * sensible_heat_flux_BlowAir


def greenhouse_air_to_thermal_screen_vapor_flux(setpoints: Setpoints, states: States):
    air_t = states.air_t
    thermal_screen_t = states.thermal_screen_t
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    air_vapor_pressure = saturation_vapor_pressure(air_t)
    thScr_vapor_pressure = saturation_vapor_pressure(thermal_screen_t)
    return air_to_obj_vapor_flux(air_vapor_pressure, thScr_vapor_pressure, HEC_AirThScr)


def greenhouse_air_to_above_thermal_screen_vapor_flux(states: States, setpoints: Setpoints, weather: Weather):
    air_t = states.air_t
    above_thermal_screen_t = states.above_thermal_screen_t
    thScr_air_flux_rate = thermal_screen_air_flux_rate(setpoints, states, weather)
    f_AirTop = thScr_air_flux_rate
    above_thermal_screen_vapor_pressure = saturation_vapor_pressure(above_thermal_screen_t)
    air_vapor_pressure = saturation_vapor_pressure(air_t)
    return general_vapor_flux(f_AirTop, air_vapor_pressure, above_thermal_screen_vapor_pressure, air_t, above_thermal_screen_t)


def greenhouse_air_to_outdoor_vapor_flux(states: States, setpoints: Setpoints, weather: Weather):
    air_t = states.air_t
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = total_side_vent_rate + f_VentForced
    air_vapor_pressure = saturation_vapor_pressure(air_t)
    outdoor_vapor_pressure = weather.outdoor_vapor_pressure
    outdoor_t = weather.outdoor_t
    return general_vapor_flux(f_AirOut, air_vapor_pressure, outdoor_vapor_pressure, air_t, outdoor_t)


def greenhouse_air_to_mechanical_cooling_vapor_flux(states: States, setpoints: Setpoints):
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    air_vapor_pressure = saturation_vapor_pressure(states.air_t)
    MechCool_vapor_pressure = saturation_vapor_pressure(states.mechcool_t)
    return air_to_obj_vapor_flux(air_vapor_pressure, MechCool_vapor_pressure, HEC_MechAir)


def above_thermal_screen_to_internal_cover_vapor_flux(states: States):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Construction.c_HECin
    cover_area = Coefficients.Construction.cover_area
    floor_area = Coefficients.Construction.floor_area
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * cover_area / floor_area
    above_thermal_screen_vapor_pressure = saturation_vapor_pressure(above_thermal_screen_t)
    cov_in_vapor_pressure = saturation_vapor_pressure(internal_cov_t)
    return air_to_obj_vapor_flux(above_thermal_screen_vapor_pressure, cov_in_vapor_pressure, HEC_TopCov_in)


def above_thermal_screen_to_outdoor_vapor_flux(states: States, setpoints: Setpoints, weather: Weather):
    above_thermal_screen_t = states.above_thermal_screen_t
    outdoor_t = weather.outdoor_t
    total_roof_vent_rate = total_roof_ventilation_rates(setpoints, states, weather)
    f_TopOut = total_roof_vent_rate
    above_thermal_screen_vapor_pressure = saturation_vapor_pressure(above_thermal_screen_t)
    outdoor_vapor_pressure = weather.outdoor_vapor_pressure
    return general_vapor_flux(f_TopOut, above_thermal_screen_vapor_pressure, outdoor_vapor_pressure, above_thermal_screen_t, outdoor_t)
