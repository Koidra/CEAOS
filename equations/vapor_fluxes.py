import math

from coefficients import Coefficients
from data_models import Setpoints, States, Weather
from equations.canopy_transpiration import canopy_transpiration_vapor_transfer_coefficient
from equations.utils import pad_and_fan_system_ventilation_flux, total_side_vents_ventilation_rates, \
    mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient, total_roof_ventilation_rates, \
    saturation_vapor_pressure
from equations.heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air, \
    thermal_screen_air_flux_rate


def air_to_obj_vapor_flux(VP_1: float, VP_2: float, HEC: float):
    # Equation 8.43
    return 6.4E-9 * HEC * max(VP_1 - VP_2, 0)


def differentiable_air_to_obj_vapor_flux(s_MV12: float, VP_1: float, VP_2: float, HEC: float):
    # Equation 8.44
    return 6.4E-9 * HEC * (VP_1 - VP_2) / (1 + math.exp(s_MV12 * (VP_1 - VP_2)))


def general_vapor_flux(f_12: float, VP_1: float, VP_2: float, object_1_t: float, object_2_t: float):
    # Equation 8.45
    M_Water = Coefficients.Outside.M_Water
    R = Coefficients.Outside.R
    return M_Water * f_12 * (VP_1 / (object_1_t + 273.15) + VP_2 / (object_2_t + 273.15)) / R


def canopy_transpiration(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """The canopy transpiration
    Equation 8.47
    :return: The canopy transpiration [kg m^-2 s^-1]
    """
    VEC_CanAir = canopy_transpiration_vapor_transfer_coefficient(states, setpoints, weather)
    VP_Can = saturation_vapor_pressure(states.can_t)
    VP_Air = saturation_vapor_pressure(states.air_t)
    return VEC_CanAir * (VP_Can - VP_Air)


def fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints: Setpoints):
    # Equation 8.64
    U_Fog = setpoints.U_Fog
    phi_Fog = Coefficients.Greenhouse.ActiveClimateControl.phi_Fog
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    return U_Fog * phi_Fog / floor_surface


# Not included in GreenLight
# def pad_and_fan_to_greenhouse_air_vapor_flux(setpoints: Setpoints):
#     # Equation 8.58
#     rho_Air = air_density()
#     f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
#     eta_Pad = Constants.Greenhouse.ActiveClimateControl.eta_Pad
#     x_Pad =
#     x_Out =
#     return rho_Air * f_Pad * (eta_Pad * (x_Pad - x_Out) + x_Out)


def heat_blower_to_greenhouse_air_vapor_flux(setpoints: Setpoints):
    # Equation 8.55
    eta_HeatVap = Coefficients.Outside.eta_HeatVap
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return eta_HeatVap * sensible_heat_flux_BlowAir


def greenhouse_air_to_outdoor_vapor_flux_by_pad_fan_system(setpoints: Setpoints, states: States):
    # Equation 8.62
    f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
    M_Water = Coefficients.Outside.M_Water
    R = Coefficients.Outside.R
    air_t = states.air_t
    VP_Air = saturation_vapor_pressure(air_t)
    return f_Pad * M_Water * (VP_Air / (air_t + 273.15)) / R


def greenhouse_air_to_thermal_screen_vapor_flux(setpoints: Setpoints, states: States):
    air_t = states.air_t
    thermal_screen_t = states.thermal_screen_t
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    VP_Air = saturation_vapor_pressure(air_t)
    VP_ThScr = saturation_vapor_pressure(thermal_screen_t)
    return air_to_obj_vapor_flux(VP_Air, VP_ThScr, HEC_AirThScr)


def greenhouse_air_to_above_thermal_screen_vapor_flux(states: States, setpoints: Setpoints, weather: Weather):
    air_t = states.air_t
    above_thermal_screen_t = states.above_thermal_screen_t
    f_ThScr = thermal_screen_air_flux_rate(setpoints, states, weather)
    f_AirTop = f_ThScr
    VP_Top = saturation_vapor_pressure(above_thermal_screen_t)
    VP_Air = saturation_vapor_pressure(air_t)
    return general_vapor_flux(f_AirTop, VP_Air, VP_Top, air_t, above_thermal_screen_t)


def greenhouse_air_to_outdoor_vapor_flux(states: States, setpoints: Setpoints, weather: Weather):
    air_t = states.air_t
    f_VentSide = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = f_VentSide + f_VentForced
    VP_Air = saturation_vapor_pressure(air_t)
    VP_Out = weather.VP_Out
    outdoor_t = weather.outdoor_t
    return general_vapor_flux(f_AirOut, VP_Air, VP_Out, air_t, outdoor_t)


def greenhouse_air_to_mechanical_cooling_vapor_flux(states: States, setpoints: Setpoints):
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    VP_Air = saturation_vapor_pressure(states.air_t)
    VP_MechCool = saturation_vapor_pressure(states.mechcool_t)
    return air_to_obj_vapor_flux(VP_Air, VP_MechCool, HEC_MechAir)


def above_thermal_screen_to_internal_cover_vapor_flux(states: States):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Greenhouse.Construction.c_HECin
    cover_surface = Coefficients.Greenhouse.Construction.cover_surface
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * cover_surface / floor_surface
    VP_Top = saturation_vapor_pressure(above_thermal_screen_t)
    VP_Cov_in = saturation_vapor_pressure(internal_cov_t)
    return air_to_obj_vapor_flux(VP_Top, VP_Cov_in, HEC_TopCov_in)


def above_thermal_screen_to_outdoor_vapor_flux(states: States, setpoints: Setpoints, weather: Weather):
    above_thermal_screen_t = states.above_thermal_screen_t
    outdoor_t = weather.outdoor_t
    f_VentRoof = total_roof_ventilation_rates(setpoints, states, weather)
    f_TopOut = f_VentRoof
    VP_Top = saturation_vapor_pressure(above_thermal_screen_t)
    VP_Out = weather.VP_Out
    return general_vapor_flux(f_TopOut, VP_Top, VP_Out, above_thermal_screen_t, outdoor_t)
