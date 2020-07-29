"""The implementation of Heat fluxes equations

Based on section 8.6
"""
import math

from configs import Constants

# 8.6.1 Global, PAR and NIR heat fluxes
from data_models import Inputs, Setpoints
from equations.utils import mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient, \
    pad_and_fan_system_ventilation_flux, total_side_vents_ventilation_rates, total_roof_ventilation_rates, \
    saturation_vapour_pressure
from equations.lumped_cover_layers import double_layer_cover_transmission_coefficient, \
    double_layer_cover_reflection_coefficient, shadingscreen_PAR_transmission_coefficient, \
    lumped_cover_conductive_heat_flux, shadingscreen_NIR_transmission_coefficient, \
    shadingscreen_NIR_reflection_coefficient, roof_thermal_screen_NIR_transmission_coefficient, \
    roof_thermal_screen_NIR_reflection_coefficient, shadingscreen_PAR_reflection_coefficient, \
    roof_thermal_screen_PAR_transmission_coefficient, roof_thermal_screen_PAR_reflection_coefficient
from equations.utils import air_density
from equations.vapour_fluxes import canopy_transpiration, fogging_system_to_greenhouse_air_latent_vapour_flux, \
    air_to_obj_vapour_flux


# TODO: inline these
def lumped_cover_virtual_NIR_transmission_coefficients(rho_CovNIR):
    # Equation 8.30
    return 1 - rho_CovNIR


def floor_virtual_NIR_transmission_coefficients():
    # Equation 8.30
    rho_FlrNIR = Constants.Greenhouse.Floor.rho_FlrNIR
    return 1 - rho_FlrNIR


def canopy_virtual_NIR_transmission_coefficient(inputs: Inputs):
    # Equation 8.31
    K_NIR = Constants.Global.K_NIR
    return math.exp(-K_NIR * inputs.LAI)


def canopy_virtual_NIR_reflection_coefficient(inputs: Inputs):
    # Equation 8.32
    rho_CanNIR = Constants.Global.rho_CanNIR
    tauhat_CanNIR = canopy_virtual_NIR_transmission_coefficient(inputs)
    return rho_CanNIR * (1 - tauhat_CanNIR)


def construction_elements_global_radiation(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.36
    I_Glob = inputs.I_Glob
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerNIR = shadingscreen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerNIR = shadingscreen_NIR_reflection_coefficient(setpoints)

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_roof_ThSrcNIR = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcNIR = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    rho_CovNIR = double_layer_cover_reflection_coefficient(tau_ShSrc_ShSrcPerNIR, rho_ShSrc_ShSrcPerNIR,
                                                           rho_roof_ThSrcNIR)
    rhoHat_CanNIR = canopy_virtual_NIR_reflection_coefficient(inputs)
    rho_FlrNIR = Constants.Greenhouse.Floor.rho_FlrNIR

    tauHat_CovNIR = lumped_cover_virtual_NIR_transmission_coefficients(rho_CovNIR)
    tauHat_CanNIR = canopy_virtual_NIR_transmission_coefficient(inputs)
    tauHat_FlrNIR = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    tau_CovCanNIR = double_layer_cover_transmission_coefficient(tauHat_CovNIR, tauHat_CanNIR, rho_CovNIR,
                                                                rhoHat_CanNIR)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    rho_CovCanNIR = double_layer_cover_reflection_coefficient(tauHat_CovNIR, rho_CovNIR,
                                                              rhoHat_CanNIR)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    tau_CovCanFlrNIR = double_layer_cover_transmission_coefficient(tau_CovCanNIR, tauHat_FlrNIR, rho_CovCanNIR,
                                                                   rho_FlrNIR)  # line 389 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover, canopy and floor
    rho_CovCanFlrNIR = double_layer_cover_reflection_coefficient(tau_CovCanNIR, rho_CovCanNIR,
                                                                 rho_FlrNIR)  # line 392 / setGlAux / GreenLight

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
    a_FlrNIR = tau_CovCanFlrNIR # page 213

    tau_CovPAR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR, rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)

    eta_GlobAir = Constants.Greenhouse.Construction.eta_GlobAir
    eta_GlobPAR = Constants.Global.eta_GlobPAR
    eta_GlobNIR = Constants.Global.eta_GlobNIR
    return eta_GlobAir * I_Glob * (tau_CovPAR * eta_GlobPAR + (a_CanNIR + a_FlrNIR) * eta_GlobNIR)


def thermal_screen_FIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.39
    U_ThScr = setpoints.U_ThScr
    # TODO: check if tau_ThScrTIR is tau_ThScrFIR
    tau_ThScrFIR = Constants.Greenhouse.Thermalscreen.tau_ThScrFIR
    return 1 - U_ThScr * (1 - tau_ThScrFIR)


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


def thermal_screen_air_flux_rate(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.41
    U_ThScr = setpoints.U_ThScr
    K_ThScr = Constants.Greenhouse.Thermalscreen.K_ThScr
    air_t = inputs.air_t
    outdoor_t = inputs.outdoor_t
    rho_Air = air_density()
    R = Constants.Global.R
    M_Air = Constants.Global.M_Air
    h_Elevation = Constants.Greenhouse.Construction.h_Elevation
    pressure = 101325 * (1 - 2.5577e-5 * h_Elevation) ** 5.25588
    rho_Top = M_Air*pressure/((inputs.above_thermal_screen_t+273.15)*R)
    rho_Out = rho_Top # = rho_Top, line 715 / setGlAux / GreenLight
    rho_mean_Air = (rho_Air + rho_Out) / 2
    g = Constants.Global.g
    return U_ThScr*K_ThScr*abs(air_t-outdoor_t)**0.66 + (1-U_ThScr)*(0.5*rho_mean_Air*(1-U_ThScr)*g*abs(rho_Air-rho_Out))**0.5/rho_mean_Air


def latent_heat_fluxes(MV) -> float:
    """The latent heat flux
    Equation 8.42
    :param float MV: the vapour flux from object 1 to object 2
    :return: The latent heat flux from object 1 to object 2 [W m^-2]
    """
    return Constants.Global.delta_H * MV


def sensible_heat_flux_between_canopy_and_air(inputs: Inputs):
    can_t = inputs.can_t
    air_t = inputs.air_t
    HEC_CanAir = 2 * Constants.Global.alpha_LeafAir * inputs.LAI
    return convective_and_conductive_heat_fluxes(HEC_CanAir, can_t, air_t)


def latent_heat_flux_between_canopy_and_air(inputs: Inputs):
    MV_CanAir = canopy_transpiration(inputs)
    return latent_heat_fluxes(MV_CanAir)


def sensible_heat_flux_between_greenhouse_air_and_outdoor_by_pad_fan_system(setpoints: Setpoints, inputs: Inputs):
    # Equation 8.61
    f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
    rho_Air = air_density()
    c_pAir = Constants.Global.c_pAir
    air_t = inputs.air_t
    return f_Pad * rho_Air * c_pAir * air_t

# Not included in GreenLight
# def sensible_heat_flux_between_pad_and_greenhouse_air(setpoints: Setpoints, inputs: Inputs):
#     # Equation 8.60
#     f_Pad = pad_and_fan_system_ventilation_flux(setpoints)
#     rho_Air = air_density()
#     c_pAir = Constants.Global.c_pAir
#     outdoor_t = inputs.outdoor_t
#     delta_H = Constants.Global.delta_H
#     eta_Pad = Constants.Greenhouse.ActiveClimateControl.eta_Pad
#     x_Pad =
#     x_Out =
#     return f_Pad * (rho_Air * c_pAir * outdoor_t - delta_H * rho_Air * (eta_Pad * (x_Pad - x_Out)))


def sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints: Setpoints, inputs: Inputs):
    air_t = inputs.air_t
    mechcool_t = inputs.mechcool_t
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, inputs)
    return convective_and_conductive_heat_fluxes(HEC_MechAir, mechcool_t, air_t)


def sensible_heat_flux_between_heating_pipe_and_greenhouse_air(inputs: Inputs):
    pipe_t = inputs.pipe_t
    air_t = inputs.air_t
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    HEC_PipeAir = 1.99 * math.pi * phi_Pipe_e * l_Pipe * abs(pipe_t - air_t) ** 0.32
    return convective_and_conductive_heat_fluxes(HEC_PipeAir, pipe_t, air_t)


def sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints: Setpoints):
    # Equation 8.53
    U_Blow = setpoints.U_Blow
    P_Blow = Constants.Greenhouse.ActiveClimateControl.P_Blow
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    return U_Blow * P_Blow / A_Flr


def sensible_heat_flux_between_buffer_and_greenhouse_air(inputs: Inputs):
    # Equation 8.57
    HEC_PasAir = Constants.Greenhouse.ActiveClimateControl.HEC_PasAir
    soil_3_t = inputs.soil_j_t[2] # third layer
    air_t = inputs.air_t
    return HEC_PasAir * (soil_3_t - air_t)

def sensible_heat_flux_between_floor_and_greenhouse_air(inputs: Inputs):
    floor_t = inputs.floor_t
    air_t = inputs.air_t
    if floor_t > air_t:
        HEC_AirFlr = 1.7 * (floor_t - air_t) ** 0.33
    else:
        HEC_AirFlr = 1.3 * (air_t - floor_t) ** 0.25
    return convective_and_conductive_heat_fluxes(HEC_AirFlr, air_t, floor_t)


def sensible_heat_flux_between_thermal_screen_and_greenhouse_air(inputs: Inputs, setpoints: Setpoints):
    air_t = inputs.air_t
    thermal_screen_t = inputs.thermal_screen_t
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_AirThScr, air_t, thermal_screen_t)


def sensible_heat_flux_between_outdoor_and_greenhouse_air(inputs: Inputs, setpoints: Setpoints):
    c_pAir = Constants.Global.c_pAir
    air_t = inputs.air_t
    outdoor_t = inputs.outdoor_t
    rho_Air = air_density()
    f_VentSide = total_side_vents_ventilation_rates(setpoints, inputs)
    f_VentForced = 0 # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    HEC_AirOut = rho_Air * c_pAir * (f_VentSide + f_VentForced)
    return convective_and_conductive_heat_fluxes(HEC_AirOut, air_t, outdoor_t)


def sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(inputs: Inputs, setpoints: Setpoints):
    c_pAir = Constants.Global.c_pAir
    air_t = inputs.air_t
    rho_Air = air_density()
    above_thermal_screen_t = inputs.above_thermal_screen_t
    f_ThScr = thermal_screen_air_flux_rate(setpoints, inputs)
    HEC_AirTop = rho_Air * c_pAir * f_ThScr
    return convective_and_conductive_heat_fluxes(HEC_AirTop, air_t, above_thermal_screen_t)


def latent_heat_flux_between_fogging_and_greenhouse_air(setpoints: Setpoints):
    MV_FogAir = fogging_system_to_greenhouse_air_latent_vapour_flux(setpoints)
    return latent_heat_fluxes(MV_FogAir)


def sensible_heat_flux_between_floor_and_first_layer_soil(inputs: Inputs):
    h_Flr = Constants.Greenhouse.Floor.h_Flr
    lambda_Flr = Constants.Greenhouse.Floor.lambda_Flr
    lambda_soil = Constants.Greenhouse.Soil.lambda_soil
    floor_t = inputs.floor_t
    h_So1 = Constants.Greenhouse.Soil.h_So[0]
    HEC_FlrSo1 = 2 / (h_Flr / lambda_Flr + h_So1 / lambda_soil)
    soil_1_t = inputs.soil_j_t[0] # first layer
    return convective_and_conductive_heat_fluxes(HEC_FlrSo1, floor_t, soil_1_t)

def latent_heat_flux_between_greenhouse_air_and_thermal_screen(inputs: Inputs, setpoints: Setpoints):
    air_t = inputs.air_t
    thermal_screen_t = inputs.thermal_screen_t
    VP_Air = saturation_vapour_pressure(air_t)
    VP_ThScr = saturation_vapour_pressure(thermal_screen_t)
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    MV_AirThScr = air_to_obj_vapour_flux(VP_Air, VP_ThScr, HEC_AirThScr)
    return latent_heat_fluxes(MV_AirThScr)


def sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(inputs: Inputs, setpoints: Setpoints):
    above_thermal_screen_t = inputs.above_thermal_screen_t
    thermal_screen_t = inputs.thermal_screen_t
    HEC_ThScrTop = 1.7 * setpoints.U_ThScr * abs(thermal_screen_t - above_thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_ThScrTop, thermal_screen_t, above_thermal_screen_t)


def sensible_heat_flux_between_above_thermal_screen_and_internal_cover(inputs: Inputs, setpoints: Setpoints):
    above_thermal_screen_t = inputs.above_thermal_screen_t
    internal_cov_t = inputs.internal_cov_t
    c_HECin = Constants.Greenhouse.Construction.c_HECin
    A_Cov = Constants.Greenhouse.Construction.A_Cov
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * A_Cov / A_Flr
    return convective_and_conductive_heat_fluxes(HEC_TopCov_in, above_thermal_screen_t, internal_cov_t)

def sensible_heat_flux_between_above_thermal_screen_and_outdoor(inputs: Inputs, setpoints: Setpoints):
    c_pAir = Constants.Global.c_pAir
    above_thermal_screen_t = inputs.above_thermal_screen_t
    outdoor_t = inputs.outdoor_t
    rho_Air = air_density()
    f_VentRoof = total_roof_ventilation_rates(setpoints, inputs)
    HEC_TopOut = rho_Air * c_pAir * f_VentRoof
    return convective_and_conductive_heat_fluxes(HEC_TopOut, above_thermal_screen_t, outdoor_t)

def latent_heat_flux_between_above_thermal_screen_and_internal_cover(inputs: Inputs):
    above_thermal_screen_t = inputs.above_thermal_screen_t
    internal_cov_t = inputs.internal_cov_t
    c_HECin = Constants.Greenhouse.Construction.c_HECin
    A_Cov = Constants.Greenhouse.Construction.A_Cov
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    VP_Top = saturation_vapour_pressure(above_thermal_screen_t)
    VP_Cov_in = saturation_vapour_pressure(internal_cov_t)
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * A_Cov / A_Flr
    MV_TopCov_in = air_to_obj_vapour_flux(VP_Top, VP_Cov_in, HEC_TopCov_in)
    return latent_heat_fluxes(MV_TopCov_in)

def sensible_heat_flux_between_internal_cover_and_external_cover(inputs: Inputs, setpoints: Setpoints):
    internal_cov_t = inputs.internal_cov_t
    external_cov_t = inputs.external_cov_t
    HEC_Cov_in_Cov_e = lumped_cover_conductive_heat_flux(setpoints) # Note: line 819 / setGlAux / GreenLight
    return convective_and_conductive_heat_fluxes(HEC_Cov_in_Cov_e, internal_cov_t, external_cov_t)


def sensible_heat_flux_between_external_cover_and_outdoor(inputs: Inputs):
    A_Cov = Constants.Greenhouse.Construction.A_Cov
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    c_HECout_1 = Constants.Greenhouse.Construction.c_HECout_1
    c_HECout_2 = Constants.Greenhouse.Construction.c_HECout_2
    external_cov_t = inputs.external_cov_t
    outdoor_t = inputs.outdoor_t
    v_Wind = inputs.v_Wind
    c_HECout_3 = Constants.Greenhouse.Construction.c_HECout_3
    HEC_Cov_e_Out = A_Cov * (c_HECout_1 + c_HECout_2 * v_Wind ** c_HECout_3) / A_Flr
    return convective_and_conductive_heat_fluxes(HEC_Cov_e_Out, external_cov_t, outdoor_t)


def heat_flux_to_heating_pipe(U, P, A):
    # Equation 8.56
    return U * P / A