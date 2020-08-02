"""
Note:
    PAR: Photosynthetically active radiation
    FIR: Far infrared radiation
    NIR: Near infrared radiation
"""
from equations.heat_fluxes import *
from equations.lumped_cover_layers import *


def canopy_PAR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
    """The PAR absorbed by the canopy
    Equation 8.26
    :return: The PAR absorbed by the canopy [W m^-2]
    """
    return canopy_PAR_absorbed_from_greenhouse_cover(states, setpoints, weather) + canopy_PAR_absorbed_from_greenhouse_floor(states, setpoints, weather)


def canopy_PAR_absorbed_from_greenhouse_cover(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.27
    radiation_flux_PARGh = PAR_above_canopy(states, setpoints, weather)
    canopy_PAR_reflection_coefficient = Coefficients.Outside.canopy_PAR_reflection_coefficient
    canopy_PAR_extinction_coefficient = Coefficients.Outside.canopy_PAR_extinction_coefficient
    return radiation_flux_PARGh * (1 - canopy_PAR_reflection_coefficient) * (1 - math.exp(-canopy_PAR_extinction_coefficient * states.leaf_area_index))


def canopy_PAR_absorbed_from_greenhouse_floor(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.29
    radiation_flux_PARGh = PAR_above_canopy(states, setpoints, weather)
    canopy_PAR_extinction_coefficient = Coefficients.Outside.canopy_PAR_extinction_coefficient
    floor_PAR_extinction_coefficient = Coefficients.Outside.floor_PAR_extinction_coefficient
    floor_PAR_reflection_coefficient = Coefficients.Greenhouse.Floor.floor_PAR_reflection_coefficient
    canopy_PAR_reflection_coefficient = Coefficients.Outside.canopy_PAR_reflection_coefficient
    return radiation_flux_PARGh * (1 - math.exp(-canopy_PAR_extinction_coefficient * states.leaf_area_index)) * floor_PAR_reflection_coefficient * (1 - canopy_PAR_reflection_coefficient) * (1 - math.exp(-floor_PAR_extinction_coefficient * states.leaf_area_index))


def floor_NIR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.34
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_NIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coefficient = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coefficient, shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)
    virtual_NIR_reflection_canopy_coef = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coefficient = Coefficients.Greenhouse.Floor.floor_NIR_reflection_coefficient

    virtual_NIR_transmission_cover_coef = lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coefficient)
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    virtual_NIR_transmission_floor_coef = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    cover_can_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(virtual_NIR_transmission_cover_coef, virtual_NIR_transmission_canopy_coef, cover_NIR_reflection_coefficient, virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    cover_can_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(virtual_NIR_transmission_cover_coef, cover_NIR_reflection_coefficient, virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    cover_can_floor_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(cover_can_NIR_transmission_coefficient, virtual_NIR_transmission_floor_coef, cover_can_NIR_reflection_coefficient, floor_NIR_reflection_coefficient)  # line 389 / setGlAux / GreenLight

    # NIR absorption coefficient of the floor
    NIR_absorption_floor_coef = cover_can_floor_NIR_transmission_coefficient  # page 213
    ratio_GlobAir = Coefficients.Greenhouse.Construction.ratio_GlobAir
    ratio_GlobNIR = Coefficients.Outside.ratio_GlobNIR
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * NIR_absorption_floor_coef * ratio_GlobNIR * outdoor_global_rad


def floor_PAR_absorbed(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.35
    floor_PAR_reflection_coefficient = Coefficients.Greenhouse.Floor.floor_PAR_reflection_coefficient
    canopy_PAR_extinction_coefficient = Coefficients.Outside.canopy_PAR_extinction_coefficient
    radiation_flux_PARGh = PAR_above_canopy(states, setpoints, weather)
    return (1 - floor_PAR_reflection_coefficient) * math.exp(-canopy_PAR_extinction_coefficient * states.leaf_area_index) * radiation_flux_PARGh


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
    shScr_PAR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coefficient = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coefficient = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coefficient, roof_thScr_PAR_transmission_coefficient, shScr_PAR_reflection_coefficient, roof_thScr_PAR_reflection_coefficient)
    ratio_GlobAir = Coefficients.Greenhouse.Construction.ratio_GlobAir
    ratio_GlobPAR = Coefficients.Outside.ratio_GlobPAR
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * cover_PAR_transmission_coefficient * ratio_GlobPAR * outdoor_global_rad


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
    shScr_NIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coefficient = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coefficient, shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)
    virtual_NIR_reflection_canopy_coef = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coefficient = Coefficients.Greenhouse.Floor.floor_NIR_reflection_coefficient

    virtual_NIR_transmission_cover_coef = lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coefficient)
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    virtual_NIR_transmission_floor_coef = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    cover_can_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(virtual_NIR_transmission_cover_coef, virtual_NIR_transmission_canopy_coef, cover_NIR_reflection_coefficient, virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover and canopy
    cover_can_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(virtual_NIR_transmission_cover_coef, cover_NIR_reflection_coefficient, virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    # NIR transmission coefficient of the cover, canopy and floor
    cover_can_floor_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(cover_can_NIR_transmission_coefficient, virtual_NIR_transmission_floor_coef, cover_can_NIR_reflection_coefficient, floor_NIR_reflection_coefficient)  # line 389 / setGlAux / GreenLight
    # NIR reflection coefficient of the cover, canopy and floor
    cover_can_floor_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(cover_can_NIR_transmission_coefficient, cover_can_NIR_reflection_coefficient, floor_NIR_reflection_coefficient)  # line 392 / setGlAux / GreenLight

    # NIR absorption coefficient of the canopy
    NIR_absorption_canopy_coef = 1 - cover_can_floor_NIR_transmission_coefficient - cover_can_floor_NIR_reflection_coefficient  # page 213
    ratio_GlobAir = Coefficients.Greenhouse.Construction.ratio_GlobAir
    ratio_GlobNIR = Coefficients.Outside.ratio_GlobNIR
    outdoor_global_rad = weather.outdoor_global_rad
    return (1 - ratio_GlobAir) * NIR_absorption_canopy_coef * ratio_GlobNIR * outdoor_global_rad


def FIR_from_pipe_to_canopy(states: States):
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coefficient = Coefficients.Greenhouse.Heating.pipe_FIR_emission_coefficient
    can_FIR_emission_coefficient = Coefficients.Outside.can_FIR_emission_coefficient
    F_PipeCan = 0.49 * (1 - math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index))
    sigma = Coefficients.Outside.sigma
    pipe_t = states.pipe_t
    can_t = states.can_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, can_FIR_emission_coefficient, F_PipeCan, sigma, pipe_t, can_t)


def FIR_from_canopy_to_internal_cover(states: States, setpoints: Setpoints):
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    sigma = Coefficients.Outside.sigma
    can_t = states.can_t
    internal_cov_t = states.internal_cov_t
    can_FIR_emission_coefficient = Coefficients.Outside.can_FIR_emission_coefficient
    A_Can = 1 - math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coefficient - cover_FIR_reflection_coefficient  # = a_CovFIR, line 271 / setGlAux
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_CanCov_in = tau_U_ThScrFIR
    return net_far_infrared_radiation_fluxes(A_Can, can_FIR_emission_coefficient, epsilon_Cov, F_CanCov_in, sigma, can_t, internal_cov_t)


def FIR_from_canopy_to_floor(states: States):
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    sigma = Coefficients.Outside.sigma
    can_t = states.can_t
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    can_FIR_emission_coefficient = Coefficients.Outside.can_FIR_emission_coefficient
    A_Can = 1 - math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    floor_FIR_emission_coefficient = Coefficients.Greenhouse.Floor.floor_FIR_emission_coefficient
    F_CanFlr = 1 - 0.49 * math.pi * pipe_length * phi_external_pipe
    floor_t = states.floor_t
    return net_far_infrared_radiation_fluxes(A_Can, can_FIR_emission_coefficient, floor_FIR_emission_coefficient, F_CanFlr, sigma, can_t, floor_t)


def FIR_from_canopy_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    sky_FIR_emission_coefficient = Coefficients.Outside.sky_FIR_emission_coefficient
    can_FIR_emission_coefficient = Coefficients.Outside.can_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    can_t = states.can_t
    A_Can = 1 - math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    F_CanSky = cover_FIR_transmission_coefficient * tau_U_ThScrFIR
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Can, can_FIR_emission_coefficient, sky_FIR_emission_coefficient, F_CanSky, sigma, can_t, sky_t)


def FIR_from_canopy_to_thermal_screen(states: States, setpoints: Setpoints):
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    can_FIR_emission_coefficient = Coefficients.Outside.can_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    thScr_FIR_emission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_FIR_emission_coefficient
    F_CanThScr = setpoints.U_ThScr
    can_t = states.can_t
    A_Can = 1 - math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    thermal_screen_t = states.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Can, can_FIR_emission_coefficient, thScr_FIR_emission_coefficient, F_CanThScr, sigma, can_t, thermal_screen_t)


def FIR_from_heating_pipe_to_floor(states: States):
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Greenhouse.Floor.floor_FIR_emission_coefficient
    pipe_FIR_emission_coefficient = Coefficients.Greenhouse.Heating.pipe_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    pipe_t = states.pipe_t
    floor_t = states.floor_t
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    F_PipeFlr = 0.49
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, floor_FIR_emission_coefficient, F_PipeFlr, sigma, pipe_t, floor_t)


def FIR_from_floor_to_internal_cover(states: States, setpoints: Setpoints):
    A_Flr = 1
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Greenhouse.Floor.floor_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coefficient - cover_FIR_reflection_coefficient  # = a_CovFIR, line 271 / setGlAux
    floor_t = states.floor_t
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_FlrCov_in = tau_U_ThScrFIR * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    internal_cov_t = states.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_Flr, floor_FIR_emission_coefficient, epsilon_Cov, F_FlrCov_in, sigma, floor_t, internal_cov_t)


def FIR_from_floor_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_Flr = 1
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    sky_FIR_emission_coefficient = Coefficients.Outside.sky_FIR_emission_coefficient
    floor_FIR_emission_coefficient = Coefficients.Greenhouse.Floor.floor_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    floor_t = states.floor_t
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    tau_CovFIR = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient,
                                                             roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight

    F_FlrSky = tau_CovFIR * tau_U_ThScrFIR * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Flr, floor_FIR_emission_coefficient, sky_FIR_emission_coefficient, F_FlrSky, sigma, floor_t, sky_t)


def FIR_from_floor_to_thermal_screen(states: States, setpoints: Setpoints):
    A_Flr = 1
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    floor_FIR_emission_coefficient = Coefficients.Greenhouse.Floor.floor_FIR_emission_coefficient
    thScr_FIR_emission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    floor_t = states.floor_t
    F_FlrThScr = setpoints.U_ThScr * (1 - 0.49 * math.pi * pipe_length * phi_external_pipe) * math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    thermal_screen_t = states.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Flr, floor_FIR_emission_coefficient, thScr_FIR_emission_coefficient, F_FlrThScr, sigma, floor_t, thermal_screen_t)


def FIR_from_heating_pipe_to_thermal_screen(states: States, setpoints: Setpoints):
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coefficient = Coefficients.Greenhouse.Heating.pipe_FIR_emission_coefficient
    thScr_FIR_emission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    F_PipeThScr = setpoints.U_ThScr * 0.49 * math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    pipe_t = states.pipe_t
    thermal_screen_t = states.thermal_screen_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, thScr_FIR_emission_coefficient, F_PipeThScr, sigma, pipe_t, thermal_screen_t)


def FIR_from_thermal_screen_to_internal_cover(states: States, setpoints: Setpoints):
    A_ThScr = 1
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coefficient - cover_FIR_reflection_coefficient  # = a_CovFIR, line 271 / setGlAux
    thScr_FIR_emission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_FIR_emission_coefficient
    sigma = Coefficients.Outside.sigma
    F_ThScrCov_in = setpoints.U_ThScr
    thermal_screen_t = states.thermal_screen_t
    internal_cov_t = states.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_ThScr, thScr_FIR_emission_coefficient, epsilon_Cov, F_ThScrCov_in, sigma, thermal_screen_t, internal_cov_t)


def FIR_from_thermal_screen_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    A_ThScr = 1
    sky_FIR_emission_coefficient = Coefficients.Outside.sky_FIR_emission_coefficient
    thScr_FIR_emission_coefficient = Coefficients.Greenhouse.Thermalscreen.thScr_FIR_emission_coefficient
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    sigma = Coefficients.Outside.sigma
    F_ThScrSky = cover_FIR_transmission_coefficient * setpoints.U_ThScr
    thermal_screen_t = states.thermal_screen_t
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_ThScr, thScr_FIR_emission_coefficient, sky_FIR_emission_coefficient, F_ThScrSky, sigma, thermal_screen_t, sky_t)


def FIR_from_heating_pipe_to_internal_cover(states: States, setpoints: Setpoints):
    sigma = Coefficients.Outside.sigma
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coefficient - cover_FIR_reflection_coefficient  # = a_CovFIR, line 271 / setGlAux
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coefficient = Coefficients.Greenhouse.Heating.pipe_FIR_emission_coefficient
    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeCov_in = tau_U_ThScrFIR * 0.49 * math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    pipe_t = states.pipe_t
    internal_cov_t = states.internal_cov_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, epsilon_Cov, F_PipeCov_in, sigma, pipe_t, internal_cov_t)


def cover_global_radiation(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.37
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_PAR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coefficient = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coefficient = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # Vanthoor PAR transmission coefficient of the lumped cover
    cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coefficient, roof_thScr_PAR_transmission_coefficient, shScr_PAR_reflection_coefficient, roof_thScr_PAR_reflection_coefficient)
    # Vanthoor PAR reflection coefficient of the lumped cover
    cover_PAR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_PAR_transmission_coefficient, shScr_PAR_reflection_coefficient, roof_thScr_PAR_reflection_coefficient)

    shScr_NIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_transmission_coefficient = roof_thermal_screen_NIR_transmission_coefficient(setpoints)
    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coefficient = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR transmission coefficient of the lumped cover
    cover_NIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_NIR_transmission_coefficient, roof_thScr_NIR_transmission_coefficient, shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)
    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coefficient, shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)

    PAR_absorption_cover_coef = absorption_coefficient(cover_PAR_transmission_coefficient, cover_PAR_reflection_coefficient)
    NIR_absorption_cover_coef = absorption_coefficient(cover_NIR_transmission_coefficient, cover_NIR_reflection_coefficient)
    ratio_GlobPAR = Coefficients.Outside.ratio_GlobPAR
    ratio_GlobNIR = Coefficients.Outside.ratio_GlobNIR
    outdoor_global_rad = weather.outdoor_global_rad
    return (PAR_absorption_cover_coef * ratio_GlobPAR + NIR_absorption_cover_coef * ratio_GlobNIR) * outdoor_global_rad


def FIR_from_external_cover_to_sky(states: States, weather: Weather):
    A_Cov_e = 1
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 260 / setGlAux / GreenLight
    epsilon_Cov = 1 - cover_FIR_transmission_coefficient - cover_FIR_reflection_coefficient  # = a_CovFIR, line 271 / setGlAux
    sky_FIR_emission_coefficient = Coefficients.Outside.sky_FIR_emission_coefficient
    F_Cov_e_Sky = 1
    sigma = Coefficients.Outside.sigma
    external_cov_t = states.external_cov_t
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Cov_e, epsilon_Cov, sky_FIR_emission_coefficient, F_Cov_e_Sky, sigma, external_cov_t, sky_t)


def FIR_from_heating_pipe_to_sky(states: States, setpoints: Setpoints, weather: Weather):
    pipe_length = Coefficients.Greenhouse.Heating.pipe_length
    phi_external_pipe = Coefficients.Greenhouse.Heating.phi_external_pipe
    canopy_FIR_extinction_coefficient = Coefficients.Outside.canopy_FIR_extinction_coefficient
    A_Pipe = math.pi * pipe_length * phi_external_pipe
    pipe_FIR_emission_coefficient = Coefficients.Greenhouse.Heating.pipe_FIR_emission_coefficient
    sky_FIR_emission_coefficient = Coefficients.Outside.sky_FIR_emission_coefficient
    shScr_FIR_transmission_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_transmission_coefficient
    shScr_FIR_reflection_coefficient = Coefficients.Greenhouse.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coefficient = Coefficients.Greenhouse.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coefficient, roof_FIR_transmission_coefficient, shScr_FIR_reflection_coefficient, roof_FIR_reflection_coefficient)  # line 255 / setGlAux / GreenLight

    tau_U_ThScrFIR = thermal_screen_FIR_transmission_coefficient(setpoints)
    F_PipeSky = cover_FIR_transmission_coefficient * tau_U_ThScrFIR * 0.49 * math.exp(-canopy_FIR_extinction_coefficient * states.leaf_area_index)
    sigma = Coefficients.Outside.sigma
    pipe_t = states.pipe_t
    sky_t = weather.sky_t
    return net_far_infrared_radiation_fluxes(A_Pipe, pipe_FIR_emission_coefficient, sky_FIR_emission_coefficient, F_PipeSky, sigma, pipe_t, sky_t)


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
