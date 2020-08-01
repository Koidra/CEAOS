import math

from coefficients import Coefficients
from data_models import States, Setpoints, Weather
from equations.lumped_cover_layers import shadingscreen_PAR_transmission_coefficient, \
    shadingscreen_PAR_reflection_coefficient, roof_thermal_screen_PAR_transmission_coefficient, \
    roof_thermal_screen_PAR_reflection_coefficient, double_layer_cover_transmission_coefficient, \
    shadingscreen_NIR_transmission_coefficient, \
    shadingscreen_NIR_reflection_coefficient, roof_thermal_screen_NIR_transmission_coefficient, \
    roof_thermal_screen_NIR_reflection_coefficient
from equations.utils import air_density, saturation_vapor_pressure


def canopy_transpiration_vapor_transfer_coefficient(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """The vapor transfer coefficient of the canopy transpiration
    Equation 8.48
    :return: The vapor transfer coefficient of the canopy transpiration [kg m^-2  Pa^-1 s^-1]
    """
    evaporation_latent_heat = Coefficients.Outside.evaporation_latent_heat
    density_air = air_density()
    c_pAir = Coefficients.Outside.c_pAir
    gamma = Coefficients.Outside.gamma
    boundary_layer_resistance = Coefficients.Outside.boundary_layer_resistance
    stomatal_resistance = canopy_stomatal_resistance(states, setpoints, weather)
    return 2 * density_air * c_pAir * states.leaf_area_index / (evaporation_latent_heat * gamma * (boundary_layer_resistance + stomatal_resistance))


def canopy_stomatal_resistance(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """
    The stomatal resistance of the canopy for vapor transport
    Equation 8.49
    :return: The stomatal resistance of the canopy [s m^-1]
    """
    return Coefficients.Outside.min_canopy_transpiration_resistance * resistance_factor(states, setpoints, weather, 'above_canopy_global_radiation') * resistance_factor(states, setpoints, weather, 'CO2_Air') * resistance_factor(states, setpoints, weather, 'VP')


def resistance_factor(states: States, setpoints: Setpoints, weather: Weather, type) -> float:
    """
    The resistance factors
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    if type == 'above_canopy_global_radiation':
        outdoor_global_rad = weather.outdoor_global_rad
        ratio_GlobAir = Coefficients.Greenhouse.Construction.ratio_GlobAir
        ratio_GlobPAR = Coefficients.Outside.ratio_GlobPAR
        ratio_GlobNIR = Coefficients.Outside.ratio_GlobNIR
        eta_LampPAR = Coefficients.Greenhouse.Lamp.eta_LampPAR
        eta_LampNIR = Coefficients.Greenhouse.Lamp.eta_LampNIR
        qLampIn = Coefficients.Greenhouse.Lamp.thetaLampMax * setpoints.U_Lamp  # line 308 / setGlAux / GreenLight
        eta_IntLampPAR = Coefficients.Greenhouse.Interlight.eta_IntLampPAR
        eta_IntLampNIR = Coefficients.Greenhouse.Interlight.eta_IntLampNIR
        qIntLampIn = Coefficients.Greenhouse.Interlight.thetaIntLampMax * setpoints.U_IntLamp
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
        cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR,
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
        rCanSun = (1 - ratio_GlobAir) * outdoor_global_rad * (ratio_GlobPAR * cover_PAR_transmission_coefficient + ratio_GlobNIR * tau_CovNIR)
        # Global radiation above the canopy from the lamps
        rCanLamp = (eta_LampPAR + eta_LampNIR) * qLampIn
        # Global radiation to the canopy from the interlight lamps
        rCanIntLamp = (eta_IntLampPAR + eta_IntLampNIR) * qIntLampIn
        # Global radiation above the canopy
        above_canopy_global_radiation = rCanSun + rCanLamp + rCanIntLamp  # Note: line 338 / setGlAux / GreenLight
        return (above_canopy_global_radiation + Coefficients.Outside.c_evap1) / (above_canopy_global_radiation + Coefficients.Outside.c_evap2)
    elif type == 'CO2_Air':
        c_evap3 = smoothed_transpiration_parameters(nth=3, setpoints=setpoints, states=states, weather=weather)
        CO2_Air = states.CO2_Air
        return 1 + c_evap3(Coefficients.Outside.eta_mg_ppm * CO2_Air - 200) ** 2
    elif type == 'VP':
        c_evap4 = smoothed_transpiration_parameters(nth=4, setpoints=setpoints, states=states, weather=weather)
        canopy_vapor_pressure = saturation_vapor_pressure(states.can_t)
        air_vapor_pressure = saturation_vapor_pressure(states.air_t)
        return 1 + c_evap4(canopy_vapor_pressure - air_vapor_pressure) ** 2


def differentiable_switch(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.51
    outdoor_global_rad = weather.outdoor_global_rad
    ratio_GlobAir = Coefficients.Greenhouse.Construction.ratio_GlobAir
    ratio_GlobPAR = Coefficients.Outside.ratio_GlobPAR
    ratio_GlobNIR = Coefficients.Outside.ratio_GlobNIR
    eta_LampPAR = Coefficients.Greenhouse.Lamp.eta_LampPAR
    eta_LampNIR = Coefficients.Greenhouse.Lamp.eta_LampNIR
    qLampIn = Coefficients.Greenhouse.Lamp.thetaLampMax * setpoints.U_Lamp  # line 308 / setGlAux / GreenLight
    eta_IntLampPAR = Coefficients.Greenhouse.Interlight.eta_IntLampPAR
    eta_IntLampNIR = Coefficients.Greenhouse.Interlight.eta_IntLampNIR
    qIntLampIn = Coefficients.Greenhouse.Interlight.thetaIntLampMax * setpoints.U_IntLamp
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
    cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(tau_ShSrc_ShSrcPerPAR, tau_roof_ThSrcPAR, rho_ShSrc_ShSrcPerPAR, rho_roof_ThSrcPAR)

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
    rCanSun = (1 - ratio_GlobAir) * outdoor_global_rad * (ratio_GlobPAR * cover_PAR_transmission_coefficient + ratio_GlobNIR * tau_CovNIR)
    # Global radiation above the canopy from the lamps
    rCanLamp = (eta_LampPAR + eta_LampNIR) * qLampIn
    # Global radiation to the canopy from the interlight lamps
    rCanIntLamp = (eta_IntLampPAR + eta_IntLampNIR) * qIntLampIn
    # Global radiation above the canopy
    above_canopy_global_radiation = rCanSun + rCanLamp + rCanIntLamp  # Note: line 338 / setGlAux / GreenLight
    return 1 / (1 + math.exp(Coefficients.Outside.s_r_s(above_canopy_global_radiation - Coefficients.Outside.rad_canopy_setpoint)))


def smoothed_transpiration_parameters(nth: int, states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.52
    S_r_s = differentiable_switch(states, setpoints, weather)
    if nth == 3:
        return Coefficients.Outside.c_night_evap3 * (1 - S_r_s) + Coefficients.Outside.c_night_evap3 * S_r_s
    else:
        return Coefficients.Outside.c_night_evap4 * (1 - S_r_s) + Coefficients.Outside.c_night_evap4 * S_r_s
