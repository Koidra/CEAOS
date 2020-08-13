"""
Note:
    PAR: Photosynthetically active radiation
    FIR: Far infrared radiation
    NIR: Near infrared radiation
"""
from climate.electrical_input import lamp_electrical_input
from climate.heat_fluxes import *
from climate.lumped_cover_layers import *
from constants import *
from data_models import States, Weather


def canopy_PAR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
    """The PAR absorbed by the canopy
    Equation 8.26
    :return: The PAR absorbed by the canopy [W m^-2]
    """
    return canopy_PAR_absorbed_from_greenhouse_cover(states, setpoints, weather) + \
           canopy_PAR_absorbed_from_greenhouse_floor(states, setpoints, weather)


def canopy_PAR_absorbed_from_greenhouse_cover(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.27
    radiation_flux_PARGh = PAR_above_canopy(states, setpoints, weather)
    return radiation_flux_PARGh * (1 - CANOPY_PAR_REFLECTION_COEF) * \
           (1 - math.exp(-CANOPY_PAR_EXTINCTION_COEF * states.leaf_area_index))


def canopy_PAR_absorbed_from_greenhouse_floor(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.29
    radiation_flux_PARGh = PAR_above_canopy(states, setpoints, weather)
    floor_PAR_reflection_coef = Coefficients.Floor.floor_PAR_reflection_coefficient
    return radiation_flux_PARGh * (1 - math.exp(-CANOPY_PAR_EXTINCTION_COEF * states.leaf_area_index)) * \
           floor_PAR_reflection_coef * (1 - CANOPY_PAR_REFLECTION_COEF) * \
           (1 - math.exp(-FLOOR_PAR_EXTINCTION_COEF * states.leaf_area_index))


def floor_NIR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.34
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_NIR_transmission_coef = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coef = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coef = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coef,
                                                                          shScr_NIR_reflection_coef,
                                                                          roof_thScr_NIR_reflection_coef)
    virtual_NIR_reflection_canopy_coef = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coef = Coefficients.Floor.floor_NIR_reflection_coefficient

    virtual_NIR_transmission_cover_coef = lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coef)
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    virtual_NIR_transmission_floor_coef = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    cover_canopy_NIR_transmission_coef = double_layer_cover_transmission_coefficient(
        virtual_NIR_transmission_cover_coef, virtual_NIR_transmission_canopy_coef, cover_NIR_reflection_coef,
        virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    cover_canopy_NIR_reflection_coef = double_layer_cover_reflection_coefficient(virtual_NIR_transmission_cover_coef,
                                                                                 cover_NIR_reflection_coef,
                                                                                 virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_transmission_coef = double_layer_cover_transmission_coefficient(
        cover_canopy_NIR_transmission_coef, virtual_NIR_transmission_floor_coef, cover_canopy_NIR_reflection_coef,
        floor_NIR_reflection_coef)  # line 389 / setGlAux / GreenLight

    # NIR absorption coefficient of the floor
    NIR_absorption_floor_coef = cover_canopy_floor_NIR_transmission_coef  # page 213
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * NIR_absorption_floor_coef * RATIO_GLOBALNIR * outdoor_global_rad


def floor_PAR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.35
    floor_PAR_reflection_coef = Coefficients.Floor.floor_PAR_reflection_coefficient
    canopy_PAR_extinction_coef = CANOPY_PAR_EXTINCTION_COEF
    radiation_flux_PARGh = PAR_above_canopy(states, setpoints, weather)
    return (1 - floor_PAR_reflection_coef) * \
           math.exp(-canopy_PAR_extinction_coef * states.leaf_area_index) * radiation_flux_PARGh


def PAR_above_canopy(states: States, setpoints: Setpoints, weather: Weather):
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
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * cover_PAR_transmission_coef * RATIO_GLOBALPAR * outdoor_global_rad


def canopy_NIR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
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
    shScr_NIR_transmission_coef = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coef = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coef = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coef,
                                                                          shScr_NIR_reflection_coef,
                                                                          roof_thScr_NIR_reflection_coef)
    virtual_NIR_reflection_canopy_coef = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coef = Coefficients.Floor.floor_NIR_reflection_coefficient

    virtual_NIR_transmission_cover_coef = lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coef)
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    virtual_NIR_transmission_floor_coef = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    cover_canopy_NIR_transmission_coef = double_layer_cover_transmission_coefficient(
        virtual_NIR_transmission_cover_coef, virtual_NIR_transmission_canopy_coef, cover_NIR_reflection_coef,
        virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    cover_canopy_NIR_reflection_coef = double_layer_cover_reflection_coefficient(virtual_NIR_transmission_cover_coef,
                                                                                 cover_NIR_reflection_coef,
                                                                                 virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_transmission_coef = double_layer_cover_transmission_coefficient(
        cover_canopy_NIR_transmission_coef, virtual_NIR_transmission_floor_coef, cover_canopy_NIR_reflection_coef,
        floor_NIR_reflection_coef)  # line 389 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_reflection_coef = double_layer_cover_reflection_coefficient(
        cover_canopy_NIR_transmission_coef, cover_canopy_NIR_reflection_coef,
        floor_NIR_reflection_coef)  # line 392 / setGlAux / GreenLight

    # NIR absorption coefficient of the canopy
    NIR_absorption_canopy_coef = 1 - cover_canopy_floor_NIR_transmission_coef - cover_canopy_floor_NIR_reflection_coef  # page 213
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * NIR_absorption_canopy_coef * RATIO_GLOBALNIR * outdoor_global_rad


def canopy_PAR_absorbed_from_lamp(states: States, setpoints: Setpoints):
    """
    Equation A19 [2]

    Returns: PAR from the lamps directly absorbed by the canopy [W m^{-2}]
    """
    return canopy_PAR_absorbed_from_lamp_Down(states, setpoints) \
         + canopy_PAR_absorbed_from_lamp_reflected_by_floor_Up(states, setpoints)


def canopy_PAR_absorbed_from_lamp_Down(states: States, setpoints: Setpoints):
    """
    Equation A17 [2]

    Returns: Total PAR from the lamps absorbed by the canopy [W m^{-2}]
    """
    radiation_flux_PARGh_Lamp = PAR_above_canopy_from_lamp(setpoints)

    return radiation_flux_PARGh_Lamp * (1 - CANOPY_PAR_REFLECTION_COEF) \
           * (1 - math.exp(-CANOPY_PAR_EXTINCTION_COEF * states.leaf_area_index))


def PAR_above_canopy_from_lamp(setpoints: Setpoints):
    """
    Equation A15 [2]

    Returns: PAR above the canopy from the lamps [W m^{-2}]
    """
    electrical_input_lamp = lamp_electrical_input(setpoints)
    return Coefficients.Lamp.lamp_electrical_input_PAR_conversion * electrical_input_lamp


def canopy_PAR_absorbed_from_lamp_reflected_by_floor_Up(states: States, setpoints: Setpoints):
    """
    Equation A18 [2]

    Returns: PAR from the lamps absorbed by the canopy after reflection from the floor [W m^{-2}]
    """
    radiation_flux_PARGh_Lamp = PAR_above_canopy_from_lamp(setpoints)
    floor_PAR_reflection_coef = Coefficients.Floor.floor_PAR_reflection_coefficient
    return radiation_flux_PARGh_Lamp * math.exp(-CANOPY_PAR_EXTINCTION_COEF * states.leaf_area_index) \
           * floor_PAR_reflection_coef * (1 - CANOPY_PAR_REFLECTION_COEF) \
           * (1 - math.exp(-FLOOR_PAR_EXTINCTION_COEF * states.leaf_area_index))


def canopy_NIR_absorbed_from_lamp(states: States, setpoints: Setpoints):
    """
    Equation A20 [2]
    Args:
        states:
        setpoints:

    Returns: NIR from the lamps absorbed by the canopy [W m^{-2}]

    """
    radiation_flux_PARGh_Lamp = PAR_above_canopy_from_lamp(setpoints)
    return Coefficients.Lamp.lamp_electrical_input_NIR_conversion * radiation_flux_PARGh_Lamp \
           * (1 - CANOPY_NIR_REFLECTION_COEF) * (1 - math.exp(-CANOPY_NIR_EXTINCTION_COEF * states.leaf_area_index))


def FIR_from_lamp_to_canopy(states: States):
    A_Lamp = Coefficients.Lamp.A_Lamp
    lamp_bottom_FIR_emission_coef = Coefficients.Lamp.bottom_lamp_emission
    F_LampCanopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Lamp, lamp_bottom_FIR_emission_coef, CANOPY_FIR_EMISSION_COEF,
                                             F_LampCanopy, states.lamp_t, states.canopy_t)


def canopy_PAR_absorbed_from_inter_lamp(setpoints: Setpoints):
    """
    Equation A24 [2]
    Args:
        setpoints:

    Returns: PAR from the interlights to the canopy lamps [W m^{-2}]
    """
    radiation_flux_PARGh_Lamp = PAR_above_canopy_from_lamp(setpoints)
    return Coefficients.Interlight.inter_lamp_electrical_input_PAR_conversion * radiation_flux_PARGh_Lamp


def canopy_NIR_absorbed_from_inter_lamp(setpoints: Setpoints):
    """
    Equation A25 [2]
    Args:
        setpoints:

    Returns: NIR from the interlight absorbed by the canopy [W m^{-2}]
    """
    radiation_flux_PARGh_Lamp = PAR_above_canopy_from_lamp(setpoints)
    return Coefficients.Interlight.inter_lamp_electrical_input_NIR_conversion * radiation_flux_PARGh_Lamp


def FIR_from_inter_lamp_to_canopy(states: States):
    A_Inter_lamp = Coefficients.Interlight.A_Inter_lamp
    inter_lamp_FIR_emission_coef = Coefficients.Interlight.inter_lamp_emission
    return net_far_infrared_radiation_fluxes(A_Inter_lamp, inter_lamp_FIR_emission_coef, CANOPY_FIR_EMISSION_COEF,
                                             1, states.intLamp_t, states.canopy_t)


def floor_PAR_absorbed_from_lamp(states: States, setpoints: Setpoints, weather: Weather):
    """
    Equation A21 [2]
    Args:
        states:
        setpoints:
        weather:

    Returns: PAR from the lamps absorbed by the floor [W m^{-2}]
    """
    radiation_flux_PARGh_Lamp = PAR_above_canopy_from_lamp(setpoints)

    return (1-Coefficients.Floor.floor_PAR_reflection_coefficient) \
           * math.exp(-CANOPY_PAR_EXTINCTION_COEF * states.leaf_area_index)\
           * radiation_flux_PARGh_Lamp


def floor_NIR_absorbed_from_lamp(states: States, setpoints: Setpoints):
    """
    Equation A22 [5]
    Args:
        states:
        setpoints:

    Returns: NIR from the lamps absorbed by the floor [W m^{-2}]
    """
    electrical_input_lamp = lamp_electrical_input(setpoints)
    return (1 - Coefficients.Floor.floor_NIR_reflection_coefficient) \
           * math.exp(-CANOPY_NIR_EXTINCTION_COEF * states.leaf_area_index) \
           * Coefficients.Lamp.lamp_electrical_input_NIR_conversion * electrical_input_lamp


def FIR_from_lamp_to_floor(states: States):
    F_LampFlr = (1 - 0.49 * math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe) \
                * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(Coefficients.Lamp.A_Lamp,
                                             Coefficients.Lamp.bottom_lamp_emission,
                                             Coefficients.Floor.floor_FIR_emission_coefficient,
                                             F_LampFlr, states.lamp_t, states.floor_t)


def FIR_from_pipe_to_canopy(states: States):
    A_Pipe = math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe
    F_PipeCanopy = 0.49 * (1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index))
    return net_far_infrared_radiation_fluxes(A_Pipe,
                                             Coefficients.Heating.pipe_FIR_emission_coefficient, CANOPY_FIR_EMISSION_COEF,
                                             F_PipeCanopy, states.pipe_t, states.canopy_t)


def FIR_from_canopy_to_internal_cover(states: States, setpoints: Setpoints):
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coef,
                                                                          shScr_FIR_reflection_coef,
                                                                          roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coef - cover_FIR_reflection_coef  # = a_CovFIR, line 271 / setGlAux
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_CanopyCov_in = tau_U_ThScrFIR
    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF, epsilon_Cov, F_CanopyCov_in,
                                             states.canopy_t, states.internal_cov_t)


def FIR_from_canopy_to_floor(states: States):
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    F_CanopyFlr = 1 - 0.49 * math.pi * pipe_length * phi_external_pipe
    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF,
                                             Coefficients.Floor.floor_FIR_emission_coefficient, F_CanopyFlr,
                                             states.canopy_t, states.floor_t)


def FIR_from_canopy_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    F_CanopySky = cover_FIR_transmission_coef * tau_U_ThScrFIR
    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF, SKY_FIR_EMISSION_COEF, F_CanopySky,
                                             states.canopy_t, weather.sky_t)


def FIR_from_canopy_to_thermal_screen(states: States, setpoints: Setpoints):
    F_CanopyThScr = setpoints.U_ThScr
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Canopy,
                                             CANOPY_FIR_EMISSION_COEF, Coefficients.Thermalscreen.thScr_FIR_emission_coefficient,
                                             F_CanopyThScr, states.canopy_t, states.thermal_screen_t)


def FIR_from_heating_pipe_to_floor(states: States):
    A_Pipe = math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe
    F_PipeFlr = 0.49
    return net_far_infrared_radiation_fluxes(A_Pipe,
                                             Coefficients.Heating.pipe_FIR_emission_coefficient,
                                             Coefficients.Floor.floor_FIR_emission_coefficient,
                                             F_PipeFlr, states.pipe_t, states.floor_t)


def FIR_from_floor_to_internal_cover(states: States, setpoints: Setpoints):
    A_Flr = 1
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    shScr_FIR_transmission_coefficient = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient,
                                                                                     roof_FIR_transmission_coefficient,
                                                                                     shScr_FIR_reflection_coefficient,
                                                                                     roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coefficient,
                                                                                 shScr_FIR_reflection_coefficient,
                                                                                 roof_FIR_reflection_coefficient)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coefficient - cover_FIR_reflection_coefficient  # = a_CovFIR, line 271 / setGlAux
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_FlrCov_in = tau_U_ThScrFIR * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) \
                  * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Flr,
                                             Coefficients.Floor.floor_FIR_emission_coefficient,
                                             epsilon_Cov,
                                             F_FlrCov_in, states.floor_t, states.internal_cov_t)


def FIR_from_floor_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_Flr = 1
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    shScr_FIR_transmission_coefficient = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Roof.roof_FIR_reflection_coefficient
    tau_CovFIR = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient,
                                                             roof_FIR_transmission_coefficient,
                                                             shScr_FIR_reflection_coefficient,
                                                             roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight

    F_FlrSky = tau_CovFIR * tau_U_ThScrFIR * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * \
               math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Flr,
                                             Coefficients.Floor.floor_FIR_emission_coefficient,
                                             SKY_FIR_EMISSION_COEF,
                                             F_FlrSky, states.floor_t, weather.sky_t)


def FIR_from_floor_to_thermal_screen(states: States, setpoints: Setpoints):
    A_Flr = 1
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    F_FlrThScr = setpoints.U_ThScr * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * \
                 math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Flr,
                                             Coefficients.Floor.floor_FIR_emission_coefficient,
                                             Coefficients.Thermalscreen.thScr_FIR_emission_coefficient, F_FlrThScr,
                                             states.floor_t, states.thermal_screen_t)


def FIR_from_heating_pipe_to_thermal_screen(states: States, setpoints: Setpoints):
    A_Pipe = math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe
    F_PipeThScr = setpoints.U_ThScr * 0.49 * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Pipe,
                                             Coefficients.Heating.pipe_FIR_emission_coefficient,
                                             Coefficients.Thermalscreen.thScr_FIR_emission_coefficient, F_PipeThScr,
                                             states.pipe_t, states.thermal_screen_t)


def FIR_from_thermal_screen_to_internal_cover(states: States, setpoints: Setpoints):
    A_ThScr = 1
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coef,
                                                                          shScr_FIR_reflection_coef,
                                                                          roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coef - cover_FIR_reflection_coef  # = a_CovFIR, line 271 / setGlAux
    F_ThScrCov_in = setpoints.U_ThScr
    return net_far_infrared_radiation_fluxes(A_ThScr,
                                             Coefficients.Thermalscreen.thScr_FIR_emission_coefficient, epsilon_Cov,
                                             F_ThScrCov_in, states.thermal_screen_t, states.internal_cov_t)


def FIR_from_thermal_screen_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_ThScr = 1
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    F_ThScrSky = cover_FIR_transmission_coef * setpoints.U_ThScr
    return net_far_infrared_radiation_fluxes(A_ThScr,
                                             Coefficients.Thermalscreen.thScr_FIR_emission_coefficient,
                                             SKY_FIR_EMISSION_COEF, F_ThScrSky,
                                             states.thermal_screen_t, weather.sky_t)


def FIR_from_heating_pipe_to_internal_cover(states: States, setpoints: Setpoints):
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coef,
                                                                          shScr_FIR_reflection_coef,
                                                                          roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coef - cover_FIR_reflection_coef  # = a_CovFIR, line 271 / setGlAux
    A_Pipe = math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeCov_in = tau_U_ThScrFIR * 0.49 * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Pipe,
                                             Coefficients.Heating.pipe_FIR_emission_coefficient, epsilon_Cov,
                                             F_PipeCov_in, states.pipe_t, states.internal_cov_t)


def cover_global_radiation(setpoints: Setpoints, weather: Weather):
    # Equation 8.37
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
    # Vanthoor PAR reflection coefficient of the lumped cover
    cover_PAR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_PAR_transmission_coef,
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
    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coef,
                                                                          shScr_NIR_reflection_coef,
                                                                          roof_thScr_NIR_reflection_coef)

    PAR_absorption_cover_coef = absorption_coefficient(cover_PAR_transmission_coef, cover_PAR_reflection_coef)
    NIR_absorption_cover_coef = absorption_coefficient(cover_NIR_transmission_coef, cover_NIR_reflection_coef)
    outdoor_global_rad = weather.outdoor_global_rad
    return (PAR_absorption_cover_coef * RATIO_GLOBALPAR + NIR_absorption_cover_coef * RATIO_GLOBALNIR) * outdoor_global_rad


def FIR_from_external_cover_to_sky(states: States, weather: Weather):
    A_Cov_e = 1
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coef,
                                                                          shScr_FIR_reflection_coef,
                                                                          roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coef - cover_FIR_reflection_coef  # = a_CovFIR, line 271 / setGlAux
    F_Cov_e_Sky = 1
    return net_far_infrared_radiation_fluxes(A_Cov_e, epsilon_Cov, SKY_FIR_EMISSION_COEF, F_Cov_e_Sky,
                                             states.external_cov_t, weather.sky_t)


def FIR_from_heating_pipe_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_Pipe = math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight

    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeSky = cover_FIR_transmission_coef * tau_U_ThScrFIR * 0.49 \
                * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Pipe,
                                             Coefficients.Heating.pipe_FIR_emission_coefficient, SKY_FIR_EMISSION_COEF,
                                             F_PipeSky, states.pipe_t, weather.sky_t)


def FIR_from_canopy_to_blackout_screen(states: States, setpoints: Setpoints):
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    F_CanopyBlScr = Coefficients.Lamp.lamp_FIR_transmission_coef * setpoints.U_BlScr

    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF,
                                             Coefficients.Blackoutscreen.blScr_FIR_emission_coef,
                                             F_CanopyBlScr, states.canopy_t, states.blScr_t)


def FIR_from_grow_pipe_to_canopy(states: States):
    A_GroPipe = math.pi * Coefficients.GrowPipe.pipe_length * Coefficients.GrowPipe.phi_external_pipe  # Surface area of pipes for floor area
    F_GroPipeCanopy = 1
    return net_far_infrared_radiation_fluxes(A_GroPipe, Coefficients.GrowPipe.groPipe_FIR_emission_coef,
                                             CANOPY_FIR_EMISSION_COEF, F_GroPipeCanopy, states.groPipe_t, states.canopy_t)


def FIR_from_floor_to_blackout_screen(states: States, setpoints: Setpoints):
    A_Flr = 1
    F_FloorBlScr = Coefficients.Lamp.lamp_FIR_transmission_coef * setpoints.U_BlScr * (1 - 0.49 * math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe) * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Flr,
                                             Coefficients.Floor.floor_FIR_emission_coefficient,
                                             Coefficients.Blackoutscreen.blScr_FIR_emission_coef,
                                             F_FloorBlScr, states.floor_t, states.blScr_t)


def FIR_from_blackout_screen_to_thermal_screen(states: States, setpoints: Setpoints):
    A_BlScr = setpoints.U_BlScr
    F_BlScrThScr = setpoints.U_ThScr
    return net_far_infrared_radiation_fluxes(A_BlScr, Coefficients.Blackoutscreen.blScr_FIR_emission_coef,
                                             Coefficients.Thermalscreen.thScr_FIR_emission_coefficient,
                                             F_BlScrThScr, states.blScr_t, states.thermal_screen_t)


def FIR_from_lamp_to_thermal_screen(states: States, setpoints: Setpoints):
    blScr_FIR_transmission_coef = 1 - setpoints.U_BlScr * (1 - Coefficients.Blackoutscreen.blScr_PAR_transmission_coef)
    F_LampThScr = setpoints.U_ThScr * blScr_FIR_transmission_coef
    return net_far_infrared_radiation_fluxes(Coefficients.Lamp.A_Lamp, Coefficients.Lamp.top_lamp_emission,
                                             Coefficients.Thermalscreen.thScr_FIR_emission_coefficient,
                                             F_LampThScr, states.lamp_t, states.thermal_screen_t)


def FIR_from_blackout_screen_to_internal_cover(states: States, setpoints: Setpoints):
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coef,
                                                                          shScr_FIR_reflection_coef,
                                                                          roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coef - cover_FIR_reflection_coef  # = a_CovFIR, line 271 / setGlAux
    F_BlScrCov_in = setpoints.U_BlScr*thermal_screen_FIR_transmission_coefficient(setpoints)
    return net_far_infrared_radiation_fluxes(1, Coefficients.Blackoutscreen.blScr_FIR_emission_coef, epsilon_Cov,
                                             F_BlScrCov_in, states.blScr_t, states.internal_cov_t)


def FIR_from_lamp_to_internal_cover(states: States, setpoints: Setpoints):
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coef,
                                                                          shScr_FIR_reflection_coef,
                                                                          roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coef - cover_FIR_reflection_coef  # = a_CovFIR, line 271 / setGlAux
    F_LampCov_in = thermal_screen_FIR_transmission_coefficient(setpoints).blackout_screen_FIR_transmission_coefficient(setpoints)
    return net_far_infrared_radiation_fluxes(Coefficients.Lamp.A_Lamp, Coefficients.Lamp.top_lamp_emission, epsilon_Cov,
                                             F_LampCov_in, states.lamp_t, states.internal_cov_t)


def FIR_from_heating_pipe_to_blackout_screen(states: States, setpoints: Setpoints):
    A_Pipe = math.pi * Coefficients.Heating.pipe_length * Coefficients.Heating.phi_external_pipe
    F_PipeBlScr = Coefficients.Lamp.lamp_FIR_transmission_coef * setpoints.U_BlScr * 0.49 \
                  * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(A_Pipe, Coefficients.Heating.pipe_FIR_emission_coefficient,
                                             Coefficients.Blackoutscreen.blScr_FIR_emission_coef,
                                             F_PipeBlScr, states.pipe_t, states.blScr_t)


def FIR_from_lamp_to_heating_pipe(states: States):
    F_LampPipe = 0.49*math.pi*Coefficients.Heating.pipe_length*Coefficients.Heating.phi_external_pipe \
                 * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    return net_far_infrared_radiation_fluxes(Coefficients.Lamp.A_Lamp,
                                             Coefficients.Lamp.bottom_lamp_emission,
                                             Coefficients.Heating.pipe_FIR_emission_coefficient,
                                             F_LampPipe, states.lamp_t, states.pipe_t)


def FIR_from_lamp_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    F_LampSky = cover_FIR_transmission_coef\
                * thermal_screen_FIR_transmission_coefficient(setpoints)\
                * blackout_screen_FIR_transmission_coefficient(setpoints)
    return net_far_infrared_radiation_fluxes(Coefficients.Lamp.A_Lamp,
                                             Coefficients.Lamp.top_lamp_emission, SKY_FIR_EMISSION_COEF,
                                             F_LampSky, states.lamp_t, weather.sky_t)


def FIR_from_blackout_screen_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    F_BlScrSky = cover_FIR_transmission_coef*setpoints.U_BlScr*tau_U_ThScrFIR
    return net_far_infrared_radiation_fluxes(1, Coefficients.Blackoutscreen.blScr_FIR_emission_coef,
                                             SKY_FIR_EMISSION_COEF,
                                             F_BlScrSky, states.blScr_t, weather.sky_t)


def FIR_from_lamp_to_blackout_screen(states: States, setpoints: Setpoints):
    F_LampBlScr = setpoints.U_BlScr
    return net_far_infrared_radiation_fluxes(Coefficients.Lamp.A_Lamp,
                                             Coefficients.Lamp.top_lamp_emission,
                                             Coefficients.Blackoutscreen.blScr_FIR_emission_coef,
                                             F_LampBlScr, states.lamp_t, states.blScr_t)


def net_far_infrared_radiation_fluxes(A_i, ep_i, ep_j, F_ij, object_i_t, object_j_t) -> float:
    """The net far infrared radiation fluxes from surface ‘i’ to ‘j’
    Equation 8.38

    :param float A_i: the surface of object ‘i’ per square meter greenhouse soil
    :param float ep_i, ep_j : the thermal infrared emission coefficients for object ‘i’ and ‘j’ respectively
    :param float F_ij: the view factor from object ‘i’ to ‘j’
    :param float object_i_t, object_j_t: the temperatures of object ‘i’ and ‘j’ respectively
    :return: The net far infrared radiation fluxes from surface ‘i’ to ‘j’ [W m^-2]
    """
    return A_i * ep_i * ep_j * F_ij * BOLTZMANN * ((object_i_t + 273.15) ** 4 - (object_j_t + 273.15) ** 4)
