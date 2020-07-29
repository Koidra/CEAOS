import math

from configs import Constants
from data_models import Inputs, Setpoints
from equations.lumped_cover_layers import shadingscreen_PAR_transmission_coefficient, \
    shadingscreen_PAR_reflection_coefficient, roof_thermal_screen_PAR_transmission_coefficient, \
    roof_thermal_screen_PAR_reflection_coefficient, double_layer_cover_transmission_coefficient, \
    double_layer_cover_reflection_coefficient, shadingscreen_NIR_transmission_coefficient, \
    shadingscreen_NIR_reflection_coefficient, roof_thermal_screen_NIR_transmission_coefficient, \
    roof_thermal_screen_NIR_reflection_coefficient
from equations.utils import air_density, saturation_vapour_pressure


def canopy_transpiration_vapour_transfer_coefficient(inputs: Inputs) -> float:
    """The vapour transfer coefficient of the canopy transpiration
    Equation 8.48
    :return: The vapour transfer coefficient of the canopy transpiration [kg m^-2  Pa^-1 s^-1]
    """
    delta_H = Constants.Global.delta_H
    rho_Air = air_density()
    c_pAir = Constants.Global.c_pAir
    gamma = Constants.Global.gamma
    r_b = Constants.Global.r_b
    r_s = canopy_stomatal_resistance(inputs)
    return 2 * rho_Air * c_pAir * inputs.LAI / (delta_H * gamma * (r_b + r_s))


def canopy_stomatal_resistance(inputs: Inputs) -> float:
    """
    The stomatal resistance of the canopy for vapour transport
    Equation 8.49
    :return: The stomatal resistance of the canopy [s m^-1]
    """
    return Constants.Global.r_s_min * resistance_factor(inputs, 'R_Can') * resistance_factor(inputs, 'CO2_Air') * resistance_factor(inputs, 'VP')


def resistance_factor(inputs: Inputs, setpoints: Setpoints, type) -> float:
    """
    The resistance factors
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    if type == 'R_Can':
        I_Glob = inputs.I_Glob
        eta_GlobAir = Constants.Greenhouse.Construction.eta_GlobAir
        eta_GlobPAR = Constants.Global.eta_GlobPAR
        eta_GlobNIR = Constants.Global.eta_GlobNIR
        eta_LampPAR = Constants.Greenhouse.Lamp.eta_LampPAR
        eta_LampNIR = Constants.Greenhouse.Lamp.eta_LampNIR
        qLampIn = Constants.Greenhouse.Lamp.thetaLampMax * setpoints.U_Lamp  # line 308 / setGlAux / GreenLight
        eta_IntLampPAR = Constants.Greenhouse.Interlight.eta_IntLampPAR
        eta_IntLampNIR = Constants.Greenhouse.Interlight.eta_IntLampNIR
        qIntLampIn = Constants.Greenhouse.Interlight.thetaIntLampMax * setpoints.U_IntLamp
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
        tau_CovPAR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR,
                                                                 rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)

        # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
        tau_ShSrc_ShSrcPerNIR = shadingscreen_NIR_transmission_coefficient(setpoints)
        # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
        rho_ShSrc_ShSrcPerNIR = shadingscreen_NIR_reflection_coefficient(setpoints)

        # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
        tau_roof_ThSrcNIR = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
        # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
        rho_roof_ThSrcNIR = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

        # Vanthoor NIR transmission coefficient of the lumped cover
        tau_CovNIR = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerNIR, tau_roof_ThSrcNIR,
                                                                 rho_ShSrc_ShSrcPerNIR, rho_roof_ThSrcNIR)

        # Global radiation above the canopy from the sun
        rCanSun = (1 - eta_GlobAir) * I_Glob * (eta_GlobPAR * tau_CovPAR + eta_GlobNIR * tau_CovNIR)
        # Global radiation above the canopy from the lamps
        rCanLamp = (eta_LampPAR + eta_LampNIR) * qLampIn
        # Global radiation to the canopy from the interlight lamps
        rCanIntLamp = (eta_IntLampPAR + eta_IntLampNIR) * qIntLampIn
        # Global radiation above the canopy
        R_Can = rCanSun + rCanLamp + rCanIntLamp  # Note: line 338 / setGlAux / GreenLight
        return (R_Can + Constants.Global.c_evap1) / (R_Can + Constants.Global.c_evap2)
    elif type == 'CO2_Air':
        c_evap3 = smoothed_transpiration_parameters(nth=3, setpoints=setpoints, inputs=inputs)
        CO2_Air = inputs.CO2_Air
        return 1 + c_evap3(Constants.Global.eta_mg_ppm * CO2_Air - 200) ** 2
    elif type == 'VP':
        c_evap4 = smoothed_transpiration_parameters(nth=4, setpoints=setpoints, inputs=inputs)
        VP_Can = saturation_vapour_pressure(inputs.can_t)
        VP_Air = saturation_vapour_pressure(inputs.air_t)
        return 1 + c_evap4(VP_Can - VP_Air) ** 2


def differentiable_switch(inputs: Inputs, setpoints: Setpoints):
    # Equation 8.51
    I_Glob = inputs.I_Glob
    eta_GlobAir = Constants.Greenhouse.Construction.eta_GlobAir
    eta_GlobPAR = Constants.Global.eta_GlobPAR
    eta_GlobNIR = Constants.Global.eta_GlobNIR
    eta_LampPAR = Constants.Greenhouse.Lamp.eta_LampPAR
    eta_LampNIR = Constants.Greenhouse.Lamp.eta_LampNIR
    qLampIn = Constants.Greenhouse.Lamp.thetaLampMax * setpoints.U_Lamp # line 308 / setGlAux / GreenLight
    eta_IntLampPAR = Constants.Greenhouse.Interlight.eta_IntLampPAR
    eta_IntLampNIR = Constants.Greenhouse.Interlight.eta_IntLampNIR
    qIntLampIn = Constants.Greenhouse.Interlight.thetaIntLampMax * setpoints.U_IntLamp
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

    # Global radiation above the canopy from the sun
    rCanSun = (1 - eta_GlobAir) * I_Glob * (eta_GlobPAR * tau_CovPAR + eta_GlobNIR * tau_CovNIR)
    # Global radiation above the canopy from the lamps
    rCanLamp = (eta_LampPAR + eta_LampNIR) * qLampIn
    # Global radiation to the canopy from the interlight lamps
    rCanIntLamp = (eta_IntLampPAR + eta_IntLampNIR) * qIntLampIn
    # Global radiation above the canopy
    R_Can = rCanSun + rCanLamp + rCanIntLamp  # Note: line 338 / setGlAux / GreenLight
    return 1 / (1 + math.exp(Constants.Global.s_r_s(R_Can - Constants.Global.R_Can_SP)))


def smoothed_transpiration_parameters(nth: int, inputs: Inputs, setpoints: Setpoints):
    # Equation 8.52
    S_r_s = differentiable_switch(inputs, setpoints)
    if nth == 3:
        return Constants.Global.c_night_evap3 * (1 - S_r_s) + Constants.Global.c_night_evap3 * S_r_s
    else:
        return Constants.Global.c_night_evap4 * (1 - S_r_s) + Constants.Global.c_night_evap4 * S_r_s
