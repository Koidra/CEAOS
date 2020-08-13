from ..constants import *
from ..data_models import ClimateStates, Weather
from .lumped_cover_layers import *
from .utils import air_density, saturation_vapor_pressure


def canopy_transpiration(states: ClimateStates, setpoints: Setpoints, weather: Weather) -> float:
    """The canopy transpiration
    Equation 8.47
    :return: The canopy transpiration [kg m^-2 s^-1]
    """
    VEC_CanopyAir = canopy_transpiration_vapor_transfer_coefficient(states, setpoints, weather)
    canopy_vp = saturation_vapor_pressure(states.t_Canopy)
    return VEC_CanopyAir * (canopy_vp - states.vapor_pressure_Air)


def canopy_transpiration_vapor_transfer_coefficient(states: ClimateStates, setpoints: Setpoints, weather: Weather) -> float:
    """The vapor transfer coefficient of the canopy transpiration
    Equation 8.48
    :return: The vapor transfer coefficient of the canopy transpiration [kg m^-2  Pa^-1 s^-1]
    """
    density_air = air_density()
    stomatal_resistance = canopy_stomatal_resistance(states, setpoints, weather)
    return 2 * density_air * C_PAIR * states.leaf_area_index / \
           (EVAPORATION_LATENT_HEAT * GAMMA * (BOUNDARY_LAYER_RESISTANCE + stomatal_resistance))


def canopy_stomatal_resistance(states: ClimateStates, setpoints: Setpoints, weather: Weather) -> float:
    """
    The stomatal resistance of the canopy for vapor transport
    Equation 8.49
    :return: The stomatal resistance of the canopy [s m^-1]
    """
    return MIN_CANOPY_TRANSPIRATION_RESISTANCE * \
           global_radiation_resistance_factor(setpoints, weather) * \
           greenhouse_co2_resistance_factor(states, setpoints, weather) * \
           vapor_pressure_resistance_factor(states, setpoints, weather)


def global_radiation_resistance_factor(setpoints: Setpoints, weather: Weather) -> float:
    """
    The resistance factors for high radiation levels
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    outdoor_global_rad = weather.outdoor_global_rad
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_PAR_transmission_coef = Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coef = Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coef = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coef = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    cover_PAR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coef,
                                                                              roof_thScr_PAR_transmission_coef,
                                                                              shScr_PAR_reflection_coef,
                                                                              roof_thScr_PAR_reflection_coef)

    shScr_NIR_transmission_coef = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coef = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_transmission_coef = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coef = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR transmission coefficient of the lumped cover
    cover_NIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_NIR_transmission_coef,
                                                                              roof_thScr_NIR_transmission_coef,
                                                                              shScr_NIR_reflection_coef,
                                                                              roof_thScr_NIR_reflection_coef)

    # Global radiation above the canopy from the sun
    rCanopySun = (1 - ratio_GlobAir) * outdoor_global_rad * \
              (RATIO_GLOBALPAR * cover_PAR_transmission_coef + RATIO_GLOBALNIR * cover_NIR_transmission_coef)
    # Global radiation above the canopy
    above_canopy_global_radiation = rCanopySun  # Note: line 338 / setGlAux / GreenLight
    return (above_canopy_global_radiation + C_EVAP1) / (above_canopy_global_radiation + C_EVAP2)


def greenhouse_co2_resistance_factor(states: ClimateStates, setpoints: Setpoints, weather: Weather) -> float:
    """
    The resistance factors for high CO2 levels
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    c_evap3 = smoothed_transpiration_parameters(nth=3, setpoints=setpoints, weather=weather)
    return 1 + c_evap3(ETA_MG_PPM * states.co2_Air - 200) ** 2


def vapor_pressure_resistance_factor(states: ClimateStates, setpoints: Setpoints, weather: Weather) -> float:
    """
    The resistance factors for large vapor pressure differences
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    c_evap4 = smoothed_transpiration_parameters(nth=4, setpoints=setpoints, weather=weather)
    canopy_vp = saturation_vapor_pressure(states.t_Canopy)
    return 1 + c_evap4(canopy_vp - states.vapor_pressure_Air) ** 2


def differentiable_switch(setpoints: Setpoints, weather: Weather):
    # Equation 8.51
    outdoor_global_rad = weather.outdoor_global_rad
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_PAR_transmission_coef = Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coef = Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coef = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coef = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    cover_PAR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coef,
                                                                              roof_thScr_PAR_transmission_coef,
                                                                              shScr_PAR_reflection_coef,
                                                                              roof_thScr_PAR_reflection_coef)

    shScr_NIR_transmission_coef = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coef = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_transmission_coef = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coef = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR transmission coefficient of the lumped cover
    cover_NIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_NIR_transmission_coef,
                                                                              roof_thScr_NIR_transmission_coef,
                                                                              shScr_NIR_reflection_coef,
                                                                              roof_thScr_NIR_reflection_coef)

    # Global radiation above the canopy from the sun
    rCanopySun = (1 - ratio_GlobAir) * outdoor_global_rad * \
                 (RATIO_GLOBALPAR * cover_PAR_transmission_coef + RATIO_GLOBALNIR * cover_NIR_transmission_coef)
    # Global radiation above the canopy
    above_canopy_global_radiation = rCanopySun  # Note: line 338 / setGlAux / GreenLight
    return 1 / (1 + math.exp(S_R_S * (above_canopy_global_radiation - RAD_CANOPY_SETPOINT)))


def smoothed_transpiration_parameters(nth: int, setpoints: Setpoints, weather: Weather):
    # Equation 8.52
    S_r_s = differentiable_switch(setpoints, weather)
    if nth == 3:
        return C_NIGHT_EVAP3 * (1 - S_r_s) + C_NIGHT_EVAP3 * S_r_s
    else:
        return C_NIGHT_EVAP4 * (1 - S_r_s) + C_NIGHT_EVAP4 * S_r_s
