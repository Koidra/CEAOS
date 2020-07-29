import math

from configs import Constants
from data_models import Setpoints, Inputs
from equations.canopy_transpiration import canopy_transpiration_vapour_transfer_coefficient
from equations.utils import pad_and_fan_system_ventilation_flux, total_side_vents_ventilation_rates, \
    mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient, total_roof_ventilation_rates, \
    saturation_vapour_pressure
from equations.heat_fluxes import sensible_heat_flux_between_direct_air_heater_and_greenhouse_air, \
    thermal_screen_air_flux_rate
from equations.utils import air_density


def air_to_obj_vapour_flux(VP_1: float, VP_2: float, HEC: float):
    # Equation 8.43
    if VP_1 < VP_2:
        return 0
    else:
        return 6.4E-9 * HEC * (VP_1 - VP_2)


def differentiable_air_to_obj_vapour_flux(s_MV12: float, VP_1: float, VP_2: float, HEC: float):
    # Equation 8.44
    return 6.4E-9 * HEC * (VP_1 - VP_2) / (1 + math.exp(s_MV12 * (VP_1 - VP_2)))


def general_vapour_flux(f_12: float, VP_1: float, VP_2: float, object_1_t: float, object_2_t: float):
    # Equation 8.45
    M_Water = Constants.Global.M_Water
    R = Constants.Global.R
    return M_Water * f_12 * (VP_1 / (object_1_t + 273.15) + VP_2 / (object_2_t + 273.15)) / R


def canopy_transpiration(inputs: Inputs) -> float:
    """The canopy transpiration
    Equation 8.47
    :return: The canopy transpiration [kg m^-2 s^-1]
    """
    VEC_CanAir = canopy_transpiration_vapour_transfer_coefficient(inputs)
    VP_Can = saturation_vapour_pressure(inputs.can_t)
    VP_Air = saturation_vapour_pressure(inputs.air_t)
    return VEC_CanAir * (VP_Can - VP_Air)


def fogging_system_to_greenhouse_air_latent_vapour_flux(setpoints: Setpoints):
    # Equation 8.64
    U_Fog = setpoints.U_Fog
    phi_Fog = Constants.Greenhouse.ActiveClimateControl.phi_Fog
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    return U_Fog * phi_Fog / A_Flr


# Not included in GreenLight
# def pad_and_fan_to_greenhouse_air_vapour_flux(setpoints: Setpoints):
#     # Equation 8.58
#     rho_Air = air_density()
#     f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
#     eta_Pad = Constants.Greenhouse.ActiveClimateControl.eta_Pad
#     x_Pad =
#     x_Out =
#     return rho_Air * f_Pad * (eta_Pad * (x_Pad - x_Out) + x_Out)


def heat_blower_to_greenhouse_air_vapour_flux(setpoints: Setpoints):
    # Equation 8.55
    eta_HeatVap = Constants.Global.eta_HeatVap
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    return eta_HeatVap * sensible_heat_flux_BlowAir


def greenhouse_air_to_outdoor_vapour_flux_by_pad_fan_system(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.62
    f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
    M_Water = Constants.Global.M_Water
    R = Constants.Global.R
    air_t = inputs.air_t
    VP_Air = saturation_vapour_pressure(air_t)
    return f_Pad * M_Water * (VP_Air / (air_t + 273.15)) / R


def greenhouse_air_to_thermal_screen_vapour_flux(setpoints: Setpoints, inputs: Inputs):
    air_t = inputs.air_t
    thermal_screen_t = inputs.thermal_screen_t
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    VP_Air = saturation_vapour_pressure(air_t)
    VP_ThScr = saturation_vapour_pressure(thermal_screen_t)
    return air_to_obj_vapour_flux(VP_Air, VP_ThScr, HEC_AirThScr)

def greenhouse_air_to_above_thermal_screen_vapour_flux(inputs: Inputs, setpoints: Setpoints):
    air_t = inputs.air_t
    above_thermal_screen_t = inputs.above_thermal_screen_t
    f_ThScr = thermal_screen_air_flux_rate(setpoints, inputs)
    f_AirTop = f_ThScr
    VP_Top = saturation_vapour_pressure(above_thermal_screen_t)
    VP_Air = saturation_vapour_pressure(air_t)
    return general_vapour_flux(f_AirTop, VP_Air, VP_Top, air_t, above_thermal_screen_t)

def greenhouse_air_to_outdoor_vapour_flux(inputs: Inputs, setpoints: Setpoints):
    air_t = inputs.air_t
    f_VentSide = total_side_vents_ventilation_rates(setpoints, inputs)
    f_VentForced = 0 # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    f_AirOut = f_VentSide + f_VentForced
    VP_Air = saturation_vapour_pressure(air_t)
    VP_Out = inputs.VP_Out
    outdoor_t = inputs.outdoor_t
    return general_vapour_flux(f_AirOut, VP_Air, VP_Out, air_t, outdoor_t)

def greenhouse_air_to_mechanical_cooling_vapour_flux(inputs: Inputs, setpoints: Setpoints):
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, inputs)
    VP_Air = saturation_vapour_pressure(inputs.air_t)
    VP_MechCool = saturation_vapour_pressure(inputs.mechcool_t)
    return air_to_obj_vapour_flux(VP_Air, VP_MechCool, HEC_MechAir)

def above_thermal_screen_to_internal_cover_vapour_flux(inputs: Inputs):
    above_thermal_screen_t = inputs.above_thermal_screen_t
    internal_cov_t = inputs.internal_cov_t
    c_HECin = Constants.Greenhouse.Construction.c_HECin
    A_Cov = Constants.Greenhouse.Construction.A_Cov
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * A_Cov / A_Flr
    VP_Top = saturation_vapour_pressure(above_thermal_screen_t)
    VP_Cov_in = saturation_vapour_pressure(internal_cov_t)
    return air_to_obj_vapour_flux(VP_Top, VP_Cov_in, HEC_TopCov_in)


def above_thermal_screen_to_outdoor_vapour_flux(inputs: Inputs, setpoints: Setpoints):
    above_thermal_screen_t = inputs.above_thermal_screen_t
    outdoor_t = inputs.outdoor_t
    f_VentRoof = total_roof_ventilation_rates(setpoints, inputs)
    f_TopOut = f_VentRoof
    VP_Top = saturation_vapour_pressure(above_thermal_screen_t)
    VP_Out = inputs.VP_Out
    return general_vapour_flux(f_TopOut, VP_Top, VP_Out, above_thermal_screen_t, outdoor_t)