"""
Note:
    PAR: Photosynthetically active radiation
    FIR: Far infrared radiation
    NIR: Near infrared radiation
"""
from climate.equations.heat_fluxes import *
from climate.equations.lumped_cover_layers import *
from constants import *


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
    cover_canopy_NIR_transmission_coef = double_layer_cover_transmission_coefficient(virtual_NIR_transmission_cover_coef,
                                                                                     virtual_NIR_transmission_canopy_coef,
                                                                                     cover_NIR_reflection_coef,
                                                                                     virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    cover_canopy_NIR_reflection_coef = double_layer_cover_reflection_coefficient(virtual_NIR_transmission_cover_coef,
                                                                                 cover_NIR_reflection_coef,
                                                                                 virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_transmission_coef = double_layer_cover_transmission_coefficient(cover_canopy_NIR_transmission_coef,
                                                                                           virtual_NIR_transmission_floor_coef,
                                                                                           cover_canopy_NIR_reflection_coef,
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
    cover_canopy_NIR_transmission_coef = double_layer_cover_transmission_coefficient(virtual_NIR_transmission_cover_coef,
                                                                                     virtual_NIR_transmission_canopy_coef,
                                                                                     cover_NIR_reflection_coef,
                                                                                     virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    cover_canopy_NIR_reflection_coef = double_layer_cover_reflection_coefficient(virtual_NIR_transmission_cover_coef,
                                                                                 cover_NIR_reflection_coef,
                                                                                 virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_transmission_coef = double_layer_cover_transmission_coefficient(cover_canopy_NIR_transmission_coef,
                                                                                           virtual_NIR_transmission_floor_coef,
                                                                                           cover_canopy_NIR_reflection_coef,
                                                                                           floor_NIR_reflection_coef)  # line 389 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_reflection_coef = double_layer_cover_reflection_coefficient(cover_canopy_NIR_transmission_coef,
                                                                                       cover_canopy_NIR_reflection_coef,
                                                                                       floor_NIR_reflection_coef)  # line 392 / setGlAux / GreenLight

    # NIR absorption coefficient of the canopy
    NIR_absorption_canopy_coef = 1 - cover_canopy_floor_NIR_transmission_coef - cover_canopy_floor_NIR_reflection_coef  # page 213
    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * NIR_absorption_canopy_coef * RATIO_GLOBALNIR * outdoor_global_rad


def FIR_from_pipe_to_canopy(states: States):
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coef = Coefficients.Heating.pipe_FIR_emission_coefficient
    canopy_FIR_emission_coef = CANOPY_FIR_EMISSION_COEF
    F_PipeCanopy = 0.49 * (1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index))
    pipe_t = states.pipe_t
    canopy_t = states.canopy_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coef, canopy_FIR_emission_coef, F_PipeCanopy,
                                             pipe_t, canopy_t)


def FIR_from_canopy_to_internal_cover(states: States, setpoints: Setpoints):
    canopy_t = states.canopy_t
    internal_cov_t = states.internal_cov_t
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
                                             canopy_t, internal_cov_t)


def FIR_from_canopy_to_floor(states: States):
    canopy_t = states.canopy_t
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    floor_FIR_emission_coefficient = Coefficients.Floor.floor_FIR_emission_coefficient
    F_CanopyFlr = 1 - 0.49 * math.pi * pipe_length * phi_external_pipe
    floor_t = states.floor_t
    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF, floor_FIR_emission_coefficient, F_CanopyFlr,
                                             canopy_t, floor_t)


def FIR_from_canopy_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    canopy_t = states.canopy_t
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
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF, SKY_FIR_EMISSION_COEF, F_CanopySky,
                                             canopy_t, sky_t)


def FIR_from_canopy_to_thermal_screen(states: States, setpoints: Setpoints):
    thScr_FIR_emission_coefficient = Coefficients.Thermalscreen.thScr_FIR_emission_coefficient
    F_CanopyThScr = setpoints.U_ThScr
    canopy_t = states.canopy_t
    A_Canopy = 1 - math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    thermal_screen_t = states.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Canopy, CANOPY_FIR_EMISSION_COEF, thScr_FIR_emission_coefficient, F_CanopyThScr,
                                             canopy_t, thermal_screen_t)


def FIR_from_heating_pipe_to_floor(states: States):
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Floor.floor_FIR_emission_coefficient
    pipe_FIR_emission_coefficient = Coefficients.Heating.pipe_FIR_emission_coefficient
    pipe_t = states.pipe_t
    floor_t = states.floor_t
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    F_PipeFlr = 0.49
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, floor_FIR_emission_coefficient, F_PipeFlr,
                                             pipe_t, floor_t)


def FIR_from_floor_to_internal_cover(states: States, setpoints: Setpoints):
    A_Flr = 1
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Floor.floor_FIR_emission_coefficient
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
    floor_t = states.floor_t
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_FlrCov_in = tau_U_ThScrFIR * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * math.exp(
        -CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    internal_cov_t = states.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_Flr, floor_FIR_emission_coefficient, epsilon_Cov, F_FlrCov_in,
                                             floor_t, internal_cov_t)


def FIR_from_floor_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_Flr = 1
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Floor.floor_FIR_emission_coefficient
    floor_t = states.floor_t
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
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Flr, floor_FIR_emission_coefficient, SKY_FIR_EMISSION_COEF, F_FlrSky,
                                             floor_t, sky_t)


def FIR_from_floor_to_thermal_screen(states: States, setpoints: Setpoints):
    A_Flr = 1
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Floor.floor_FIR_emission_coefficient
    thScr_FIR_emission_coefficient = Coefficients.Thermalscreen.thScr_FIR_emission_coefficient
    floor_t = states.floor_t
    F_FlrThScr = setpoints.U_ThScr * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * \
                 math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    thermal_screen_t = states.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Flr, floor_FIR_emission_coefficient, thScr_FIR_emission_coefficient, F_FlrThScr,
                                             floor_t, thermal_screen_t)


def FIR_from_heating_pipe_to_thermal_screen(states: States, setpoints: Setpoints):
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coefficient = Coefficients.Heating.pipe_FIR_emission_coefficient
    thScr_FIR_emission_coefficient = Coefficients.Thermalscreen.thScr_FIR_emission_coefficient
    F_PipeThScr = setpoints.U_ThScr * 0.49 * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    pipe_t = states.pipe_t
    thermal_screen_t = states.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, thScr_FIR_emission_coefficient, F_PipeThScr,
                                             pipe_t, thermal_screen_t)


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
    thScr_FIR_emission_coef = Coefficients.Thermalscreen.thScr_FIR_emission_coefficient
    F_ThScrCov_in = setpoints.U_ThScr
    thermal_screen_t = states.thermal_screen_t
    internal_cov_t = states.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_ThScr, thScr_FIR_emission_coef, epsilon_Cov, F_ThScrCov_in,
                                             thermal_screen_t, internal_cov_t)


def FIR_from_thermal_screen_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_ThScr = 1
    thScr_FIR_emission_coef = Coefficients.Thermalscreen.thScr_FIR_emission_coefficient
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    F_ThScrSky = cover_FIR_transmission_coef * setpoints.U_ThScr
    thermal_screen_t = states.thermal_screen_t
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_ThScr, thScr_FIR_emission_coef, SKY_FIR_EMISSION_COEF, F_ThScrSky,
                                             thermal_screen_t, sky_t)


def FIR_from_heating_pipe_to_internal_cover(states: States, setpoints: Setpoints):
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
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
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coef = Coefficients.Heating.pipe_FIR_emission_coefficient
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeCov_in = tau_U_ThScrFIR * 0.49 * math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    pipe_t = states.pipe_t
    internal_cov_t = states.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coef, epsilon_Cov, F_PipeCov_in,
                                             pipe_t, internal_cov_t)


def cover_global_radiation(states: States, setpoints: Setpoints, weather: Weather):
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
    external_cov_t = states.external_cov_t
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Cov_e, epsilon_Cov, SKY_FIR_EMISSION_COEF, F_Cov_e_Sky,
                                             external_cov_t, sky_t)


def FIR_from_heating_pipe_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coef = Coefficients.Heating.pipe_FIR_emission_coefficient
    shScr_FIR_transmission_coef = Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight

    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeSky = cover_FIR_transmission_coef * tau_U_ThScrFIR * 0.49 * \
                math.exp(-CANOPY_FIR_EXTINCTION_COEF * states.leaf_area_index)
    pipe_t = states.pipe_t
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coef, SKY_FIR_EMISSION_COEF, F_PipeSky,
                                             pipe_t, sky_t)


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
