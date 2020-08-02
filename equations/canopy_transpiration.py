from coefficients import Constants
from data_models import States, Weather
from equations.lumped_cover_layers import *
from equations.utils import air_density, saturation_vapor_pressure


def canopy_transpiration(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """The canopy transpiration
    Equation 8.47
    :return: The canopy transpiration [kg m^-2 s^-1]
    """
    VEC_CanAir = canopy_transpiration_vapor_transfer_coefficient(states, setpoints, weather)
    canopy_vapor_pressure = saturation_vapor_pressure(states.can_t)
    air_vapor_pressure = saturation_vapor_pressure(states.air_t)
    return VEC_CanAir * (canopy_vapor_pressure - air_vapor_pressure)


def canopy_transpiration_vapor_transfer_coefficient(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """The vapor transfer coefficient of the canopy transpiration
    Equation 8.48
    :return: The vapor transfer coefficient of the canopy transpiration [kg m^-2  Pa^-1 s^-1]
    """
    evaporation_latent_heat = Constants.evaporation_latent_heat
    density_air = air_density()
    c_pAir = Constants.c_pAir
    gamma = Constants.gamma
    boundary_layer_resistance = Constants.boundary_layer_resistance
    stomatal_resistance = canopy_stomatal_resistance(states, setpoints, weather)
    return 2 * density_air * c_pAir * states.leaf_area_index / (evaporation_latent_heat * gamma * (boundary_layer_resistance + stomatal_resistance))


def canopy_stomatal_resistance(states: States, setpoints: Setpoints, weather: Weather) -> float:
    """
    The stomatal resistance of the canopy for vapor transport
    Equation 8.49
    :return: The stomatal resistance of the canopy [s m^-1]
    """
    return Constants.min_canopy_transpiration_resistance * resistance_factor(states, setpoints, weather, 'above_canopy_global_radiation') * resistance_factor(states, setpoints, weather, 'CO2_Air') * resistance_factor(states, setpoints, weather, 'VP')


def resistance_factor(states: States, setpoints: Setpoints, weather: Weather, type) -> float:
    """
    The resistance factors
    Equation 8.50
    :return: The resistance factors [W m^-2]
    """
    if type == 'above_canopy_global_radiation':
        outdoor_global_rad = weather.outdoor_global_rad
        ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
        ratio_GlobPAR = Constants.ratio_GlobPAR
        ratio_GlobNIR = Constants.ratio_GlobNIR
        # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
        shScr_PAR_transmission_coefficient = Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
        shScr_PAR_reflection_coefficient = Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

        # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
        roof_thScr_PAR_transmission_coefficient = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
        # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
        roof_thScr_PAR_reflection_coefficient = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

        # Vanthoor PAR transmission coefficient of the lumped cover
        cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coefficient, roof_thScr_PAR_transmission_coefficient,
                                                                                         shScr_PAR_reflection_coefficient, roof_thScr_PAR_reflection_coefficient)

        shScr_NIR_transmission_coefficient = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
        shScr_NIR_reflection_coefficient = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

        # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
        roof_thScr_NIR_transmission_coefficient = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
        # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
        roof_thScr_NIR_reflection_coefficient = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

        # Vanthoor NIR transmission coefficient of the lumped cover
        cover_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_NIR_transmission_coefficient, roof_thScr_NIR_transmission_coefficient,
                                                                                         shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)

        # Global radiation above the canopy from the sun
        rCanSun = (1 - ratio_GlobAir) * outdoor_global_rad * (ratio_GlobPAR * cover_PAR_transmission_coefficient + ratio_GlobNIR * cover_NIR_transmission_coefficient)
        # Global radiation above the canopy
        above_canopy_global_radiation = rCanSun  # Note: line 338 / setGlAux / GreenLight
        return (above_canopy_global_radiation + Constants.c_evap1) / (above_canopy_global_radiation + Constants.c_evap2)
    elif type == 'CO2_Air':
        c_evap3 = smoothed_transpiration_parameters(nth=3, setpoints=setpoints, states=states, weather=weather)
        CO2_Air = states.CO2_Air
        return 1 + c_evap3(Constants.eta_mg_ppm * CO2_Air - 200) ** 2
    elif type == 'VP':
        c_evap4 = smoothed_transpiration_parameters(nth=4, setpoints=setpoints, states=states, weather=weather)
        canopy_vapor_pressure = saturation_vapor_pressure(states.can_t)
        air_vapor_pressure = saturation_vapor_pressure(states.air_t)
        return 1 + c_evap4(canopy_vapor_pressure - air_vapor_pressure) ** 2


def differentiable_switch(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.51
    outdoor_global_rad = weather.outdoor_global_rad
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    ratio_GlobPAR = Constants.ratio_GlobPAR
    ratio_GlobNIR = Constants.ratio_GlobNIR
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_PAR_transmission_coefficient = Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coefficient = Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coefficient = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coefficient = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coefficient, roof_thScr_PAR_transmission_coefficient, shScr_PAR_reflection_coefficient, roof_thScr_PAR_reflection_coefficient)

    shScr_NIR_transmission_coefficient = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_transmission_coefficient = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coefficient = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR transmission coefficient of the lumped cover
    cover_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_NIR_transmission_coefficient, roof_thScr_NIR_transmission_coefficient, shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)

    # Global radiation above the canopy from the sun
    rCanSun = (1 - ratio_GlobAir) * outdoor_global_rad * (ratio_GlobPAR * cover_PAR_transmission_coefficient + ratio_GlobNIR * cover_NIR_transmission_coefficient)
    # Global radiation above the canopy
    above_canopy_global_radiation = rCanSun  # Note: line 338 / setGlAux / GreenLight
    return 1 / (1 + math.exp(Constants.s_r_s(above_canopy_global_radiation - Constants.rad_canopy_setpoint)))


def smoothed_transpiration_parameters(nth: int, states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.52
    S_r_s = differentiable_switch(states, setpoints, weather)
    if nth == 3:
        return Constants.c_night_evap3 * (1 - S_r_s) + Constants.c_night_evap3 * S_r_s
    else:
        return Constants.c_night_evap4 * (1 - S_r_s) + Constants.c_night_evap4 * S_r_s
