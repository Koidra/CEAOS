from data_models import States, Weather
from climate.equations.lumped_cover_layers import *
from climate.equations.utils import air_density, saturation_vapor_pressure
from constants import *


def canopy_transpiration(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """The canopy transpiration
    Equation 8.47
    :return: The canopy transpiration [kg m^-2 s^-1]
    """
    VEC_CanopyAir = canopy_transpiration_vapor_transfer_coefficient(states, setpoints, weather)
    canopy_vp = saturation_vapor_pressure(states.canopy_t)
    air_vp = saturation_vapor_pressure(states.air_t)
    return VEC_CanopyAir * (canopy_vp - air_vp)


def canopy_transpiration_vapor_transfer_coefficient(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """The vapor transfer coefficient of the canopy transpiration
    Equation 8.48
    :return: The vapor transfer coefficient of the canopy transpiration [kg m^-2  Pa^-1 s^-1]
    """
    density_air = air_density()
    stomatal_resistance = canopy_stomatal_resistance(states, setpoints, weather)
    return 2 * density_air * C_PAIR * states.leaf_area_index / \
           (EVAPORATION_LATENT_HEAT * GAMMA * (BOUNDARY_LAYER_RESISTANCE + stomatal_resistance))


def canopy_stomatal_resistance(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """
    The stomatal resistance of the canopy for vapor transport
    Equation 8.49
    :return: The stomatal resistance of the canopy [s m^-1]
    """
    return MIN_CANOPY_TRANSPIRATION_RESISTANCE * \
           resistance_factor(states, setpoints, weather, 'above_canopy_global_radiation') * \
           resistance_factor(states, setpoints, weather, 'CO2_air') * \
           resistance_factor(states, setpoints, weather, 'VP')


def resistance_factor(states: States, setpoints: Setpoints, weather: Weather, type) -> float:
    """
    The resistance factors
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    if type == 'above_canopy_global_radiation':
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
    elif type == 'CO2_air':
        c_evap3 = smoothed_transpiration_parameters(nth=3, setpoints=setpoints, states=states, weather=weather)
        CO2_air = states.CO2_air
        return 1 + c_evap3(ETA_MG_PPM * CO2_air - 200) ** 2
    elif type == 'VP':
        c_evap4 = smoothed_transpiration_parameters(nth=4, setpoints=setpoints, states=states, weather=weather)
        canopy_vp = saturation_vapor_pressure(states.canopy_t)
        air_vp = saturation_vapor_pressure(states.air_t)
        return 1 + c_evap4(canopy_vp - air_vp) ** 2


def differentiable_switch(states: States, setpoints: Setpoints, weather: Weather):
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


def smoothed_transpiration_parameters(nth: int, states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.52
    S_r_s = differentiable_switch(states, setpoints, weather)
    if nth == 3:
        return C_NIGHT_EVAP3 * (1 - S_r_s) + C_NIGHT_EVAP3 * S_r_s
    else:
        return C_NIGHT_EVAP4 * (1 - S_r_s) + C_NIGHT_EVAP4 * S_r_s
