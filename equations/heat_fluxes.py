"""The implementation of Heat fluxes equations

Based on section 8.6
"""
import math

from coefficients import Coefficients
# 8.6.1 Global, PAR and NIR heat fluxes
from data_models import States, Setpoints, Weather
from equations.lumped_cover_layers import double_layer_cover_transmission_coefficient, \
    double_layer_cover_reflection_coefficient, shadingscreen_PAR_transmission_coefficient, \
    lumped_cover_conductive_heat_flux, shadingscreen_NIR_transmission_coefficient, \
    shadingscreen_NIR_reflection_coefficient, roof_thermal_screen_NIR_reflection_coefficient, \
    shadingscreen_PAR_reflection_coefficient, \
    roof_thermal_screen_PAR_transmission_coefficient, roof_thermal_screen_PAR_reflection_coefficient
from equations.utils import air_density
from equations.utils import mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient, \
    pad_and_fan_system_ventilation_flux, total_side_vents_ventilation_rates, total_roof_ventilation_rates, \
    saturation_vapor_pressure
from equations.vapor_fluxes import canopy_transpiration, fogging_system_to_greenhouse_air_latent_vapor_flux, \
    air_to_obj_vapor_flux


# TODO: inline these
def lumped_cover_virtual_NIR_transmission_coefficients(rho_CovNIR):
    # Equation 8.30
    return 1 - rho_CovNIR


def floor_virtual_NIR_transmission_coefficients():
    # Equation 8.30
    floor_NIR_reflection_coefficient = Coefficients.Greenhouse.Floor.floor_NIR_reflection_coefficient
    return 1 - floor_NIR_reflection_coefficient


def canopy_virtual_NIR_transmission_coefficient(states: States):
    # Equation 8.31
    canopy_NIR_extinction_coefficient = Coefficients.Outside.canopy_NIR_extinction_coefficient
    return math.exp(-canopy_NIR_extinction_coefficient * states.LAI)


def canopy_virtual_NIR_reflection_coefficient(states: States):
    # Equation 8.32
    canopy_NIR_reflection_coefficient = Coefficients.Outside.canopy_NIR_reflection_coefficient
    tauhat_CanNIR = canopy_virtual_NIR_transmission_coefficient(states)
    return canopy_NIR_reflection_coefficient * (1 - tauhat_CanNIR)


def construction_elements_global_radiation(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.36
    I_Glob = weather.I_Glob
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerNIR = shadingscreen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerNIR = shadingscreen_NIR_reflection_coefficient(setpoints)

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcNIR = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    rho_CovNIR = double_layer_cover_reflection_coefficient(tau_ShSrc_ShSrcPerNIR, rho_ShSrc_ShSrcPerNIR, rho_roof_ThSrcNIR)
    rhoHat_CanNIR = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coefficient = Coefficients.Greenhouse.Floor.floor_NIR_reflection_coefficient

    tauHat_CovNIR = lumped_cover_virtual_NIR_transmission_coefficients(rho_CovNIR)
    tauHat_CanNIR = canopy_virtual_NIR_transmission_coefficient(states)
    tauHat_FlrNIR = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    tau_CovCanNIR = double_layer_cover_transmission_coefficient(tauHat_CovNIR, tauHat_CanNIR, rho_CovNIR, rhoHat_CanNIR)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    rho_CovCanNIR = double_layer_cover_reflection_coefficient(tauHat_CovNIR, rho_CovNIR, rhoHat_CanNIR)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    tau_CovCanFlrNIR = double_layer_cover_transmission_coefficient(tau_CovCanNIR, tauHat_FlrNIR, rho_CovCanNIR, floor_NIR_reflection_coefficient)  # line 389 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover, canopy and floor
    rho_CovCanFlrNIR = double_layer_cover_reflection_coefficient(tau_CovCanNIR, rho_CovCanNIR, floor_NIR_reflection_coefficient)  # line 392 / setGlAux / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerPAR = shadingscreen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerPAR = shadingscreen_PAR_reflection_coefficient(setpoints)

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_roof_ThSrcPAR = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcPAR = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # NIR absorption coefficient of the canopy
    a_CanNIR = 1 - tau_CovCanFlrNIR - rho_CovCanFlrNIR # page 213
    # NIR absorption coefficient of the floor
    floor_surfaceNIR = tau_CovCanFlrNIR  # page 213

    tau_CovPAR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR, rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)

    eta_GlobAir = Coefficients.Greenhouse.Construction.eta_GlobAir
    ratio_GlobPAR = Coefficients.Outside.ratio_GlobPAR
    ratio_GlobNIR = Coefficients.Outside.ratio_GlobNIR
    return eta_GlobAir * I_Glob * (tau_CovPAR * ratio_GlobPAR + (a_CanNIR + floor_surfaceNIR) * ratio_GlobNIR)


def thermal_screen_FIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.39
    U_ThScr = setpoints.U_ThScr
    thScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_FIR_transmission_coefficient
    return 1 - U_ThScr * (1 - thScr_FIR_transmission_coefficient)


# 8.6.3 Convection and conduction
def convective_and_conductive_heat_fluxes(HEC, T1, T2) -> float:
    """Convective and conductive heat fluxes
    Equation 8.40

    :param float HEC: the heat exchange coefficient between object 1 and 2
    :param float T1: the temperature of object 1
    :param float T2: the temperature of object 2
    :return: The heat flow from object 1 to object 2 [W m^-2]
    """
    return HEC * (T1 - T2)


def thermal_screen_air_flux_rate(setpoints: Setpoints, states: States, weather: Weather):
    # Equation 8.41
    U_ThScr = setpoints.U_ThScr
    thScr_flux_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_flux_coefficient
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    density_air = air_density()
    R = Coefficients.Outside.R
    M_Air = Coefficients.Outside.M_Air
    elevation_height = Coefficients.Greenhouse.Construction.elevation_height
    pressure = 101325 * (1 - 2.5577e-5 * elevation_height) ** 5.25588
    rho_Top = M_Air*pressure/((states.above_thermal_screen_t+273.15)*R)
    rho_Out = rho_Top # = rho_Top, line 715 / setGlAux / GreenLight
    rho_mean_Air = (density_air + rho_Out) / 2
    g = Coefficients.Outside.g
    return U_ThScr * thScr_flux_coefficient * abs(air_t-outdoor_t) ** 0.66 + (1-U_ThScr) * (0.5 * rho_mean_Air * (1-U_ThScr) * g * abs(density_air - rho_Out)) ** 0.5 / rho_mean_Air


def latent_heat_fluxes(MV) -> float:
    """The latent heat flux
    Equation 8.42
    :param float MV: the vapor flux from object 1 to object 2
    :return: The latent heat flux from object 1 to object 2 [W m^-2]
    """
    return Coefficients.Outside.evaporation_latent_heat * MV


def sensible_heat_flux_between_canopy_and_air(states: States):
    can_t = states.can_t
    air_t = states.air_t
    HEC_CanAir = 2 * Coefficients.Outside.canopy_air_convective_heat_exchange_coefficient * states.LAI
    return convective_and_conductive_heat_fluxes(HEC_CanAir, can_t, air_t)


def latent_heat_flux_between_canopy_and_air(states: States, setpoints: Setpoints, weather: Weather):
    MV_CanAir = canopy_transpiration(states, setpoints, weather)
    return latent_heat_fluxes(MV_CanAir)


def sensible_heat_flux_between_greenhouse_air_and_outdoor_by_pad_fan_system(setpoints: Setpoints, states: States):
    # Equation 8.61
    f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
    density_air = air_density()
    c_pAir = Coefficients.Outside.c_pAir
    air_t = states.air_t
    return f_Pad * density_air * c_pAir * air_t

# Not included in GreenLight
# def sensible_heat_flux_between_pad_and_greenhouse_air(setpoints: Setpoints, states: States):
#     # Equation 8.60
#     f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
#     density_air = air_density()
#     c_pAir = Constants.Global.c_pAir
#     outdoor_t = states.outdoor_t
#     evaporation_latent_heat = Constants.Global.evaporation_latent_heat
#     eta_Pad = Constants.Greenhouse.ActiveClimateControl.eta_Pad
#     x_Pad =
#     x_Out =
#     return f_Pad * (density_air * c_pAir * outdoor_t - evaporation_latent_heat * density_air * (eta_Pad * (x_Pad - x_Out)))


def sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints: Setpoints, states: States):
    air_t = states.air_t
    mechcool_t = states.mechcool_t
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    return convective_and_conductive_heat_fluxes(HEC_MechAir, mechcool_t, air_t)


def sensible_heat_flux_between_heating_pipe_and_greenhouse_air(states: States):
    pipe_t = states.pipe_t
    air_t = states.air_t
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    HEC_PipeAir = 1.99 * math.pi * phi_external_pipe * pipe_length * abs(pipe_t - air_t) ** 0.32
    return convective_and_conductive_heat_fluxes(HEC_PipeAir, pipe_t, air_t)


def sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints: Setpoints):
    # Equation 8.53
    U_Blow = setpoints.U_Blow
    P_Blow = Coefficients.Greenhouse.ActiveClimateControl.P_Blow
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    return U_Blow * P_Blow / floor_surface


def sensible_heat_flux_between_buffer_and_greenhouse_air(states: States):
    # Equation 8.57
    HEC_PasAir = Coefficients.Greenhouse.ActiveClimateControl.HEC_PasAir
    soil_3_t = states.soil_j_t[2] # third layer
    air_t = states.air_t
    return HEC_PasAir * (soil_3_t - air_t)

def sensible_heat_flux_between_floor_and_greenhouse_air(states: States):
    floor_t = states.floor_t
    air_t = states.air_t
    if floor_t > air_t:
        HEC_AirFlr = 1.7 * (floor_t - air_t) ** 0.33
    else:
        HEC_AirFlr = 1.3 * (air_t - floor_t) ** 0.25
    return convective_and_conductive_heat_fluxes(HEC_AirFlr, air_t, floor_t)


def sensible_heat_flux_between_thermal_screen_and_greenhouse_air(states: States, setpoints: Setpoints):
    air_t = states.air_t
    thermal_screen_t = states.thermal_screen_t
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_AirThScr, air_t, thermal_screen_t)


def sensible_heat_flux_between_outdoor_and_greenhouse_air(states: States, setpoints: Setpoints, weather: Weather):
    c_pAir = Coefficients.Outside.c_pAir
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    density_air = air_density()
    f_VentSide = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0 # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    HEC_AirOut = density_air * c_pAir * (f_VentSide + f_VentForced)
    return convective_and_conductive_heat_fluxes(HEC_AirOut, air_t, outdoor_t)


def sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(states: States, setpoints: Setpoints, weather: Weather):
    c_pAir = Coefficients.Outside.c_pAir
    air_t = states.air_t
    density_air = air_density()
    above_thermal_screen_t = states.above_thermal_screen_t
    f_ThScr = thermal_screen_air_flux_rate(setpoints, states, weather)
    HEC_AirTop = density_air * c_pAir * f_ThScr
    return convective_and_conductive_heat_fluxes(HEC_AirTop, air_t, above_thermal_screen_t)


def latent_heat_flux_between_fogging_and_greenhouse_air(setpoints: Setpoints):
    MV_FogAir = fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints)
    return latent_heat_fluxes(MV_FogAir)


def sensible_heat_flux_between_floor_and_first_layer_soil(states: States):
    floor_thickness = Coefficients.Greenhouse.Floor.floor_thickness
    floor_heat_conductivity = Coefficients.Greenhouse.Floor.floor_heat_conductivity
    soil_heat_conductivity = Coefficients.Greenhouse.Soil.soil_heat_conductivity
    floor_t = states.floor_t
    h_So1 = Coefficients.Greenhouse.Soil.soil_thicknesses[0]
    HEC_FlrSo1 = 2 / (floor_thickness / floor_heat_conductivity + h_So1 / soil_heat_conductivity)
    soil_1_t = states.soil_j_t[0]  # first layer
    return convective_and_conductive_heat_fluxes(HEC_FlrSo1, floor_t, soil_1_t)


def latent_heat_flux_between_greenhouse_air_and_thermal_screen(states: States, setpoints: Setpoints):
    air_t = states.air_t
    thermal_screen_t = states.thermal_screen_t
    VP_Air = saturation_vapor_pressure(air_t)
    VP_ThScr = saturation_vapor_pressure(thermal_screen_t)
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    MV_AirThScr = air_to_obj_vapor_flux(VP_Air, VP_ThScr, HEC_AirThScr)
    return latent_heat_fluxes(MV_AirThScr)


def sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(states: States, setpoints: Setpoints):
    above_thermal_screen_t = states.above_thermal_screen_t
    thermal_screen_t = states.thermal_screen_t
    HEC_ThScrTop = 1.7 * setpoints.U_ThScr * abs(thermal_screen_t - above_thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_ThScrTop, thermal_screen_t, above_thermal_screen_t)


def sensible_heat_flux_between_above_thermal_screen_and_internal_cover(states: States, setpoints: Setpoints):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Greenhouse.Construction.c_HECin
    cover_surface = Coefficients.Greenhouse.Construction.cover_surface
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * cover_surface / floor_surface
    return convective_and_conductive_heat_fluxes(HEC_TopCov_in, above_thermal_screen_t, internal_cov_t)


def sensible_heat_flux_between_above_thermal_screen_and_outdoor(states: States, setpoints: Setpoints, weather: Weather):
    c_pAir = Coefficients.Outside.c_pAir
    above_thermal_screen_t = states.above_thermal_screen_t
    outdoor_t = weather.outdoor_t
    density_air = air_density()
    f_VentRoof = total_roof_ventilation_rates(setpoints, states, weather)
    HEC_TopOut = density_air * c_pAir * f_VentRoof
    return convective_and_conductive_heat_fluxes(HEC_TopOut, above_thermal_screen_t, outdoor_t)


def latent_heat_flux_between_above_thermal_screen_and_internal_cover(states: States):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Greenhouse.Construction.c_HECin
    cover_surface = Coefficients.Greenhouse.Construction.cover_surface
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    VP_Top = saturation_vapor_pressure(above_thermal_screen_t)
    VP_Cov_in = saturation_vapor_pressure(internal_cov_t)
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * cover_surface / floor_surface
    MV_TopCov_in = air_to_obj_vapor_flux(VP_Top, VP_Cov_in, HEC_TopCov_in)
    return latent_heat_fluxes(MV_TopCov_in)

def sensible_heat_flux_between_internal_cover_and_external_cover(states: States, setpoints: Setpoints):
    internal_cov_t = states.internal_cov_t
    external_cov_t = states.external_cov_t
    HEC_Cov_in_Cov_e = lumped_cover_conductive_heat_flux(setpoints) # Note: line 819 / setGlAux / GreenLight
    return convective_and_conductive_heat_fluxes(HEC_Cov_in_Cov_e, internal_cov_t, external_cov_t)


def sensible_heat_flux_between_external_cover_and_outdoor(states: States, weather: Weather):
    cover_surface = Coefficients.Greenhouse.Construction.cover_surface
    floor_surface = Coefficients.Greenhouse.Construction.floor_surface
    c_HECout_1 = Coefficients.Greenhouse.Construction.c_HECout_1
    c_HECout_2 = Coefficients.Greenhouse.Construction.c_HECout_2
    external_cov_t = states.external_cov_t
    outdoor_t = weather.outdoor_t
    v_Wind = weather.v_Wind
    c_HECout_3 = Coefficients.Greenhouse.Construction.c_HECout_3
    HEC_Cov_e_Out = cover_surface * (c_HECout_1 + c_HECout_2 * v_Wind ** c_HECout_3) / floor_surface
    return convective_and_conductive_heat_fluxes(HEC_Cov_e_Out, external_cov_t, outdoor_t)


def heat_flux_to_heating_pipe(U, P, A):
    # Equation 8.56
    return U * P / A