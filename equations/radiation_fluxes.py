"""
Note:
    PAR: Photosynthetically active radiation
    FIR: Far infrared radiation
    NIR: Near infrared radiation
"""
import math

from configs import Constants
from data_models import Inputs, Setpoints
from equations.heat_fluxes import thermal_screen_FIR_transmission_coefficient, \
    lumped_cover_virtual_NIR_transmission_coefficients, canopy_virtual_NIR_transmission_coefficient, \
    floor_virtual_NIR_transmission_coefficients, canopy_virtual_NIR_reflection_coefficient
from equations.lumped_cover_layers import double_layer_cover_transmission_coefficient, \
    shadingscreen_PAR_transmission_coefficient, shadingscreen_PAR_reflection_coefficient, \
    roof_thermal_screen_PAR_transmission_coefficient, roof_thermal_screen_PAR_reflection_coefficient, \
    double_layer_cover_reflection_coefficient, shadingscreen_NIR_transmission_coefficient, \
    shadingscreen_NIR_reflection_coefficient, roof_thermal_screen_NIR_transmission_coefficient, \
    roof_thermal_screen_NIR_reflection_coefficient, absorption_coefficient


def canopy_PAR_absorbed(inputs: Inputs, setpoints: Setpoints):
    """The PAR absorbed by the canopy
    Equation 8.26
    :return: The PAR absorbed by the canopy [W m^-2]
    """
    return canopy_PAR_absorbed_from_greenhouse_cover(inputs, setpoints) + canopy_PAR_absorbed_from_greenhouse_floor(inputs, setpoints)


def canopy_PAR_absorbed_from_greenhouse_cover(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.27
    radiation_flux_PARGh = PAR_above_canopy(inputs, setpoints)
    rho_CanPAR = Constants.Global.rho_CanPAR
    K_1PAR = Constants.Global.K_1PAR
    return radiation_flux_PARGh * (1 - rho_CanPAR) * (1 - math.exp(-K_1PAR * inputs.LAI))


def canopy_PAR_absorbed_from_greenhouse_floor(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.29
    radiation_flux_PARGh = PAR_above_canopy(inputs, setpoints)
    K_1PAR = Constants.Global.K_1PAR
    K_2PAR = Constants.Global.K_2PAR
    rho_FlrPAR = Constants.Greenhouse.Floor.rho_FlrPAR
    rho_CanPAR = Constants.Global.rho_CanPAR
    return radiation_flux_PARGh * (1 - math.exp(-K_1PAR * inputs.LAI)) * rho_FlrPAR * (1 - rho_CanPAR) * (1 - math.exp(-K_2PAR * inputs.LAI))


def floor_NIR_absorbed(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.34
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerNIR = shadingscreen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerNIR = shadingscreen_NIR_reflection_coefficient(setpoints)

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcNIR = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    rho_CovNIR = double_layer_cover_reflection_coefficient(tau_ShSrc_ShSrcPerNIR, rho_ShSrc_ShSrcPerNIR, rho_roof_ThSrcNIR)
    rhoHat_CanNIR = canopy_virtual_NIR_reflection_coefficient(inputs)
    rho_FlrNIR = Constants.Greenhouse.Floor.rho_FlrNIR

    tauHat_CovNIR = lumped_cover_virtual_NIR_transmission_coefficients(rho_CovNIR)
    tauHat_CanNIR = canopy_virtual_NIR_transmission_coefficient(inputs)
    tauHat_FlrNIR = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    tau_CovCanNIR = double_layer_cover_transmission_coefficient(tauHat_CovNIR, tauHat_CanNIR, rho_CovNIR, rhoHat_CanNIR)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    rho_CovCanNIR = double_layer_cover_reflection_coefficient(tauHat_CovNIR, rho_CovNIR, rhoHat_CanNIR)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    tau_CovCanFlrNIR = double_layer_cover_transmission_coefficient(tau_CovCanNIR, tauHat_FlrNIR, rho_CovCanNIR, rho_FlrNIR)  # line 389 / setGlAux / GreenLight

    # NIR absorption coefficient of the floor
    a_FlrNIR = tau_CovCanFlrNIR # page 213
    eta_GlobAir = Constants.Greenhouse.Construction.eta_GlobAir
    eta_GlobNIR = Constants.Global.eta_GlobNIR
    I_Glob = inputs.I_Glob
    return (1 - eta_GlobAir) * a_FlrNIR * eta_GlobNIR * I_Glob


def floor_PAR_absorbed(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.35
    rho_FlrPAR = Constants.Greenhouse.Floor.rho_FlrPAR
    K_1PAR = Constants.Global.K_1PAR
    radiation_flux_PARGh = PAR_above_canopy(inputs, setpoints)
    return (1 - rho_FlrPAR) * math.exp(-K_1PAR * inputs.LAI) * radiation_flux_PARGh


def PAR_above_canopy(inputs: Inputs, setpoints: Setpoints):
    """The PAR above the canopy

    The model contains four cover layers:
    + A movable outdoor shading screen
    + A semi-permanent shading screen
    + The greenhouse roof
    + A movable indoor thermal screen

    Equation 8.28
    :return: the PAR above the canopy [W m^-2]
    """
    # TODO: need to re-verify the order of four cover layers
    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerPAR = shadingscreen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerPAR = shadingscreen_PAR_reflection_coefficient(setpoints)

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_roof_ThSrcPAR = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcPAR = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    tau_CovPAR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR, rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)
    eta_GlobAir = Constants.Greenhouse.Construction.eta_GlobAir
    eta_GlobPAR = Constants.Global.eta_GlobPAR
    I_Glob = inputs.I_Glob
    return (1 - eta_GlobAir) * tau_CovPAR * eta_GlobPAR * I_Glob


def canopy_NIR_absorbed(inputs: Inputs, setpoints: Setpoints):
    """The NIR absorbed by the canopy

    The model contains four cover layers:
    + A movable outdoor shading screen
    + A semi-permanent shading screen
    + The greenhouse roof
    + A movable indoor thermal screen

    Equation 8.33
    :return: The NIR absorbed by the canopy [W m^-2]
    """
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
    rho_CovNIR = double_layer_cover_reflection_coefficient(tau_ShSrc_ShSrcPerNIR, rho_ShSrc_ShSrcPerNIR, rho_roof_ThSrcNIR)
    rhoHat_CanNIR = canopy_virtual_NIR_reflection_coefficient(inputs)
    rho_FlrNIR = Constants.Greenhouse.Floor.rho_FlrNIR

    tauHat_CovNIR = lumped_cover_virtual_NIR_transmission_coefficients(rho_CovNIR)
    tauHat_CanNIR = canopy_virtual_NIR_transmission_coefficient(inputs)
    tauHat_FlrNIR = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    tau_CovCanNIR = double_layer_cover_transmission_coefficient(tauHat_CovNIR, tauHat_CanNIR, rho_CovNIR, rhoHat_CanNIR) # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    rho_CovCanNIR = double_layer_cover_reflection_coefficient(tauHat_CovNIR, rho_CovNIR, rhoHat_CanNIR) # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    tau_CovCanFlrNIR = double_layer_cover_transmission_coefficient(tau_CovCanNIR, tauHat_FlrNIR, rho_CovCanNIR, rho_FlrNIR) # line 389 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover, canopy and floor
    rho_CovCanFlrNIR = double_layer_cover_reflection_coefficient(tau_CovCanNIR, rho_CovCanNIR, rho_FlrNIR) # line 392 / setGlAux / GreenLight

    # NIR absorption coefficient of the canopy
    a_CanNIR = 1 - tau_CovCanFlrNIR - rho_CovCanFlrNIR # page 213
    eta_GlobAir = Constants.Greenhouse.Construction.eta_GlobAir
    eta_GlobNIR = Constants.Global.eta_GlobNIR
    I_Glob = inputs.I_Glob
    return (1 - eta_GlobAir) * a_CanNIR * eta_GlobNIR * I_Glob


def FIR_from_pipe_to_canopy(inputs: Inputs):
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    K_FIR = Constants.Global.K_FIR
    A_Pipe = math.pi * l_Pipe * phi_Pipe_e
    epsilon_Pipe = Constants.Greenhouse.Heating.epsilon_Pipe
    epsilon_Can = Constants.Global.epsilon_Can
    F_PipeCan = 0.49 * (1 - math.exp(-K_FIR * inputs.LAI))
    sigma = Constants.Global.sigma
    pipe_t = inputs.pipe_t
    can_t = inputs.can_t
    return net_far_infrared_radiation_fluxes(A_Pipe, epsilon_Pipe, epsilon_Can, F_PipeCan, sigma, pipe_t, can_t)


def FIR_from_canopy_to_internal_cover(inputs: Inputs, setpoints: Setpoints):
    K_FIR = Constants.Global.K_FIR
    sigma = Constants.Global.sigma
    can_t = inputs.can_t
    internal_cov_t = inputs.internal_cov_t
    epsilon_Can = Constants.Global.epsilon_Can
    A_Can = 1 - math.exp(-K_FIR * inputs.LAI)
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR, rho_ShScrPerFIR) # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR, rho_ShScrPerFIR)# line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR =Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR, rho_RfFIR) # line 255 / setGlAux / GreenLight
    rho_CovFIR = double_layer_cover_reflection_coefficient(tau_ShScr_ShScrPerFIR, rho_ShScr_ShScrPerFIR, rho_RfFIR) # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - tau_CovFIR - rho_CovFIR # = a_CovFIR, line 271 / setGlAux
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_CanCov_in = tau_U_ThScrFIR
    return net_far_infrared_radiation_fluxes(A_Can, epsilon_Can, epsilon_Cov, F_CanCov_in, sigma, can_t, internal_cov_t)


def FIR_from_canopy_to_floor(inputs: Inputs):
    K_FIR = Constants.Global.K_FIR
    sigma = Constants.Global.sigma
    can_t = inputs.can_t
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    epsilon_Can = Constants.Global.epsilon_Can
    A_Can = 1 - math.exp(-K_FIR * inputs.LAI)
    epsilon_Flr = Constants.Greenhouse.Floor.epsilon_FlrFIR
    F_CanFlr = 1 - 0.49 * math.pi * l_Pipe * phi_Pipe_e
    floor_t = inputs.floor_t
    return net_far_infrared_radiation_fluxes(A_Can, epsilon_Can, epsilon_Flr, F_CanFlr, sigma, can_t, floor_t)


def FIR_from_canopy_to_sky(inputs: Inputs, setpoints: Setpoints):
    K_FIR = Constants.Global.K_FIR
    epsilon_Sky = Constants.Global.epsilon_Sky
    epsilon_Can = Constants.Global.epsilon_Can
    sigma = Constants.Global.sigma
    can_t = inputs.can_t
    A_Can = 1 - math.exp(-K_FIR * inputs.LAI)
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR, rho_ShScrPerFIR) # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR, rho_ShScrPerFIR)# line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR =Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR, rho_RfFIR) # line 255 / setGlAux / GreenLight
    F_CanSky = tau_CovFIR * tau_U_ThScrFIR
    sky_t = inputs.sky_t
    return net_far_infrared_radiation_fluxes(A_Can, epsilon_Can, epsilon_Sky, F_CanSky, sigma, can_t, sky_t)


def FIR_from_canopy_to_thermal_screen(inputs: Inputs, setpoints: Setpoints):
    K_FIR = Constants.Global.K_FIR
    epsilon_Can = Constants.Global.epsilon_Can
    sigma = Constants.Global.sigma
    epsilon_ThScr = Constants.Greenhouse.Thermalscreen.epsilon_ThScrFIR
    F_CanThScr = setpoints.U_ThScr
    can_t = inputs.can_t
    A_Can = 1 - math.exp(-K_FIR * inputs.LAI)
    thermal_screen_t = inputs.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Can, epsilon_Can, epsilon_ThScr, F_CanThScr, sigma, can_t, thermal_screen_t)


def FIR_from_heating_pipe_to_floor(inputs: Inputs):
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    epsilon_Flr = Constants.Greenhouse.Floor.epsilon_FlrFIR
    epsilon_Pipe = Constants.Greenhouse.Heating.epsilon_Pipe
    sigma = Constants.Global.sigma
    pipe_t = inputs.pipe_t
    floor_t = inputs.floor_t
    A_Pipe = math.pi * l_Pipe * phi_Pipe_e
    F_PipeFlr = 0.49
    return net_far_infrared_radiation_fluxes(A_Pipe, epsilon_Pipe, epsilon_Flr, F_PipeFlr, sigma, pipe_t, floor_t)


def FIR_from_floor_to_internal_cover(inputs: Inputs, setpoints: Setpoints):
    A_Flr = 1
    K_FIR = Constants.Global.K_FIR
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    epsilon_Flr = Constants.Greenhouse.Floor.epsilon_FlrFIR
    sigma = Constants.Global.sigma
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR, rho_ShScrPerFIR) # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR, rho_ShScrPerFIR)# line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR =Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR, rho_RfFIR) # line 255 / setGlAux / GreenLight
    rho_CovFIR = double_layer_cover_reflection_coefficient(tau_ShScr_ShScrPerFIR, rho_ShScr_ShScrPerFIR, rho_RfFIR) # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - tau_CovFIR - rho_CovFIR # = a_CovFIR, line 271 / setGlAux
    floor_t = inputs.floor_t
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_FlrCov_in = tau_U_ThScrFIR * (1 - 0.49 * math.pi * l_Pipe * phi_Pipe_e) * math.exp(-K_FIR * inputs.LAI)
    internal_cov_t = inputs.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_Flr, epsilon_Flr, epsilon_Cov, F_FlrCov_in, sigma, floor_t, internal_cov_t)


def FIR_from_floor_to_sky(inputs: Inputs, setpoints: Setpoints):
    A_Flr = 1
    K_FIR = Constants.Global.K_FIR
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    epsilon_Sky = Constants.Global.epsilon_Sky
    epsilon_Flr = Constants.Greenhouse.Floor.epsilon_FlrFIR
    sigma = Constants.Global.sigma
    floor_t = inputs.floor_t
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR,
                                                                        rho_ShScrPerFIR)  # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR,
                                                                      rho_ShScrPerFIR)  # line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR = Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR,
                                                             rho_RfFIR)  # line 255 / setGlAux / GreenLight

    F_FlrSky = tau_CovFIR * tau_U_ThScrFIR * (1 - 0.49 * math.pi * l_Pipe * phi_Pipe_e) * math.exp(-K_FIR * inputs.LAI)
    sky_t = inputs.sky_t
    return net_far_infrared_radiation_fluxes(A_Flr, epsilon_Flr, epsilon_Sky, F_FlrSky, sigma, floor_t, sky_t)


def FIR_from_floor_to_thermal_screen(inputs: Inputs, setpoints: Setpoints):
    A_Flr = 1
    K_FIR = Constants.Global.K_FIR
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    epsilon_Flr = Constants.Greenhouse.Floor.epsilon_FlrFIR
    epsilon_ThScr = Constants.Greenhouse.Thermalscreen.epsilon_ThScrFIR
    sigma = Constants.Global.sigma
    floor_t = inputs.floor_t
    F_FlrThScr = setpoints.U_ThScr * (1 - 0.49 * math.pi * l_Pipe * phi_Pipe_e) * math.exp(-K_FIR * inputs.LAI)
    thermal_screen_t = inputs.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Flr, epsilon_Flr, epsilon_ThScr, F_FlrThScr, sigma, floor_t, thermal_screen_t)


def FIR_from_heating_pipe_to_thermal_screen(inputs: Inputs, setpoints: Setpoints):
    K_FIR = Constants.Global.K_FIR
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    A_Pipe = math.pi * l_Pipe * phi_Pipe_e
    epsilon_Pipe = Constants.Greenhouse.Heating.epsilon_Pipe
    epsilon_ThScr = Constants.Greenhouse.Thermalscreen.epsilon_ThScrFIR
    sigma = Constants.Global.sigma
    F_PipeThScr = setpoints.U_ThScr * 0.49 * math.exp(-K_FIR * inputs.LAI)
    pipe_t = inputs.pipe_t
    thermal_screen_t = inputs.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Pipe, epsilon_Pipe, epsilon_ThScr, F_PipeThScr, sigma, pipe_t, thermal_screen_t)


def FIR_from_thermal_screen_to_internal_cover(inputs: Inputs, setpoints: Setpoints):
    A_ThScr = 1
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR,
                                                                        rho_ShScrPerFIR)  # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR,
                                                                      rho_ShScrPerFIR)  # line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR = Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR,
                                                             rho_RfFIR)  # line 255 / setGlAux / GreenLight
    rho_CovFIR = double_layer_cover_reflection_coefficient(tau_ShScr_ShScrPerFIR, rho_ShScr_ShScrPerFIR,
                                                           rho_RfFIR)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - tau_CovFIR - rho_CovFIR  # = a_CovFIR, line 271 / setGlAux
    epsilon_ThScr = Constants.Greenhouse.Thermalscreen.epsilon_ThScrFIR
    sigma = Constants.Global.sigma
    F_ThScrCov_in = setpoints.U_ThScr
    thermal_screen_t = inputs.thermal_screen_t
    internal_cov_t = inputs.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_ThScr, epsilon_ThScr, epsilon_Cov, F_ThScrCov_in, sigma, thermal_screen_t, internal_cov_t)


def FIR_from_thermal_screen_to_sky(inputs: Inputs, setpoints: Setpoints):
    A_ThScr = 1
    epsilon_Sky = Constants.Global.epsilon_Sky
    epsilon_ThScr = Constants.Greenhouse.Thermalscreen.epsilon_ThScrFIR
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR,
                                                                        rho_ShScrPerFIR)  # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR,
                                                                      rho_ShScrPerFIR)  # line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR = Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR,
                                                             rho_RfFIR)  # line 255 / setGlAux / GreenLight
    sigma = Constants.Global.sigma
    F_ThScrSky = tau_CovFIR * setpoints.U_ThScr
    thermal_screen_t = inputs.thermal_screen_t
    sky_t = inputs.sky_t
    return net_far_infrared_radiation_fluxes(A_ThScr, epsilon_ThScr, epsilon_Sky, F_ThScrSky, sigma, thermal_screen_t, sky_t)


def FIR_from_heating_pipe_to_internal_cover(inputs: Inputs, setpoints: Setpoints):
    sigma = Constants.Global.sigma
    K_FIR = Constants.Global.K_FIR
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR,
                                                                        rho_ShScrPerFIR)  # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR,
                                                                      rho_ShScrPerFIR)  # line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR = Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR,
                                                             rho_RfFIR)  # line 255 / setGlAux / GreenLight
    rho_CovFIR = double_layer_cover_reflection_coefficient(tau_ShScr_ShScrPerFIR, rho_ShScr_ShScrPerFIR,
                                                           rho_RfFIR)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - tau_CovFIR - rho_CovFIR  # = a_CovFIR, line 271 / setGlAux
    A_Pipe = math.pi * l_Pipe * phi_Pipe_e
    epsilon_Pipe = Constants.Greenhouse.Heating.epsilon_Pipe
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeCov_in = tau_U_ThScrFIR * 0.49 * math.exp(-K_FIR * inputs.LAI)
    pipe_t = inputs.pipe_t
    internal_cov_t = inputs.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_Pipe, epsilon_Pipe, epsilon_Cov, F_PipeCov_in, sigma, pipe_t, internal_cov_t)


def cover_global_radiation(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.37
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerPAR = shadingscreen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerPAR = shadingscreen_PAR_reflection_coefficient(setpoints)

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_roof_ThSrcPAR = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcPAR = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    tau_CovPAR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR, rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)
    # Vanthoor PAR reflection coefficient of the lumped cover
    rho_CovPAR = double_layer_cover_reflection_coefficient(tau_ShSrc_ShSrcPerPAR, rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_ShSrc_ShSrcPerNIR = shadingscreen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_ShSrc_ShSrcPerNIR = shadingscreen_NIR_reflection_coefficient(setpoints)

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    tau_roof_ThSrcNIR = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    rho_roof_ThSrcNIR = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR transmission coefficient of the lumped cover
    tau_CovNIR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerNIR, tau_roof_ThSrcNIR, rho_ShSrc_ShSrcPerNIR, rho_roof_ThSrcNIR)
    # Vanthoor NIR reflection coefficient of the lumped cover
    rho_CovNIR = double_layer_cover_reflection_coefficient(tau_ShSrc_ShSrcPerNIR, rho_ShSrc_ShSrcPerNIR, rho_roof_ThSrcNIR)

    a_GhPAR = absorption_coefficient(tau_CovPAR, rho_CovPAR)
    a_GhNIR =absorption_coefficient(tau_CovNIR, rho_CovNIR)
    eta_GlobPAR = Constants.Global.eta_GlobPAR
    eta_GlobNIR = Constants.Global.eta_GlobNIR
    I_Glob = inputs.I_Glob
    return (a_GhPAR * eta_GlobPAR + a_GhNIR * eta_GlobNIR) * I_Glob


def FIR_from_external_cover_to_sky(inputs: Inputs):
    A_Cov_e = 1
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR,
                                                                        rho_ShScrPerFIR)  # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR,
                                                                      rho_ShScrPerFIR)  # line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR = Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR,
                                                             rho_RfFIR)  # line 255 / setGlAux / GreenLight
    rho_CovFIR = double_layer_cover_reflection_coefficient(tau_ShScr_ShScrPerFIR, rho_ShScr_ShScrPerFIR,
                                                           rho_RfFIR)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - tau_CovFIR - rho_CovFIR  # = a_CovFIR, line 271 / setGlAux
    epsilon_Sky = Constants.Global.epsilon_Sky
    F_Cov_e_Sky = 1
    sigma = Constants.Global.sigma
    external_cov_t = inputs.external_cov_t
    sky_t = inputs.sky_t
    return net_far_infrared_radiation_fluxes(A_Cov_e, epsilon_Cov, epsilon_Sky, F_Cov_e_Sky, sigma, external_cov_t, sky_t)


def FIR_from_heating_pipe_to_sky(inputs: Inputs, setpoints: Setpoints):
    l_Pipe = Constants.Greenhouse.Heating.l_Pipe
    phi_Pipe_e = Constants.Greenhouse.Heating.phi_Pipe_e
    K_FIR = Constants.Global.K_FIR
    A_Pipe = math.pi * l_Pipe * phi_Pipe_e
    epsilon_Pipe = Constants.Greenhouse.Heating.epsilon_Pipe
    epsilon_Sky = Constants.Global.epsilon_Sky
    tau_ShScrFIR = Constants.Greenhouse.Shadowscreen.tau_ShScrFIR
    rho_ShScrFIR = Constants.Greenhouse.Shadowscreen.rho_ShScrFIR
    tau_ShScrPerFIR = Constants.Greenhouse.Whitewash.tau_ShScrPerFIR
    rho_ShScrPerFIR = Constants.Greenhouse.Whitewash.rho_ShScrPerFIR
    tau_ShScr_ShScrPerFIR = double_layer_cover_transmission_coefficient(tau_ShScrFIR, tau_ShScrPerFIR, rho_ShScrFIR,
                                                                        rho_ShScrPerFIR)  # line 103 / setGlAux / GreenLight
    rho_ShScr_ShScrPerFIR = double_layer_cover_reflection_coefficient(tau_ShScrFIR, rho_ShScrFIR,
                                                                      rho_ShScrPerFIR)  # line 107, 111 / setGlAux / GreenLight
    tau_RfFIR = Constants.Greenhouse.Roof.tau_RfFIR
    rho_RfFIR = Constants.Greenhouse.Roof.rho_RfFIR
    tau_CovFIR = double_layer_cover_transmission_coefficient(tau_ShScr_ShScrPerFIR, tau_RfFIR, rho_ShScr_ShScrPerFIR,
                                                             rho_RfFIR)  # line 255 / setGlAux / GreenLight

    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeSky = tau_CovFIR * tau_U_ThScrFIR * 0.49 * math.exp(-K_FIR * inputs.LAI)
    sigma = Constants.Global.sigma
    pipe_t = inputs.pipe_t
    sky_t = inputs.sky_t
    return net_far_infrared_radiation_fluxes(A_Pipe, epsilon_Pipe, epsilon_Sky, F_PipeSky, sigma, pipe_t, sky_t)

def net_far_infrared_radiation_fluxes(A_i, ep_i, ep_j, F_ij, sigma, object_i_t, object_j_t) -> float:
    """The net far infrared radiation fluxes from surface ‘i’ to ‘j’
    Equation 8.38

    :param float A_i: the surface of object ‘i’ per square meter greenhouse soil
    :param float ep_i, ep_j : the thermal infrared emission coefficients for object ‘i’ and ‘j’ respectively
    :param float F_ij: the view factor from object ‘i’ to ‘j’
    :param float sigma: the Stefan Boltzmann constant
    :param float T_i, T_j: the temperatures of object ‘i’ and ‘j’ respectively
    :return: The net far infrared radiation fluxes from surface ‘i’ to ‘j’ [W m^-2]
    """
    return A_i * ep_i * ep_j * F_ij * sigma * ((object_i_t + 273.15) ** 4 - (object_j_t + 273.15) ** 4)
