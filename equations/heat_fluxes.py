"""The implementation of Heat fluxes equations

Based on section 8.6
"""

# 8.6.1 Global, PAR and NIR heat fluxes
from coefficients import Constants
from equations.canopy_transpiration import canopy_transpiration
from equations.lumped_cover_layers import *
from equations.utils import *
from equations.vapor_fluxes import fogging_system_to_greenhouse_air_latent_vapor_flux, air_to_obj_vapor_flux


# TODO: inline these
def lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coefficient):
    # Equation 8.30
    return 1 - cover_NIR_reflection_coefficient


def floor_virtual_NIR_transmission_coefficients():
    # Equation 8.30
    floor_NIR_reflection_coefficient = Coefficients.Floor.floor_NIR_reflection_coefficient
    return 1 - floor_NIR_reflection_coefficient


def canopy_virtual_NIR_transmission_coefficient(states: States):
    # Equation 8.31
    canopy_NIR_extinction_coefficient = Constants.canopy_NIR_extinction_coefficient
    return math.exp(-canopy_NIR_extinction_coefficient * states.leaf_area_index)


def canopy_virtual_NIR_reflection_coefficient(states: States):
    # Equation 8.32
    canopy_NIR_reflection_coefficient = Constants.canopy_NIR_reflection_coefficient
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    return canopy_NIR_reflection_coefficient * (1 - virtual_NIR_transmission_canopy_coef)


def construction_elements_global_radiation(states: States, setpoints: Setpoints, weather: Weather):
    # Equation 8.36
    outdoor_global_rad = weather.outdoor_global_rad
    # TODO: need to re-verify the order of four cover layers and cover-canopy-floor
    shScr_NIR_transmission_coefficient = Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coefficient = Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight

    # NIR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_NIR_reflection_coefficient = roof_thermal_screen_NIR_reflection_coefficient(setpoints)

    # Vanthoor NIR reflection coefficient of the lumped cover
    cover_NIR_reflection_coefficient = double_layer_cover_reflection_coefficient(shScr_NIR_transmission_coefficient, shScr_NIR_reflection_coefficient, roof_thScr_NIR_reflection_coefficient)
    virtual_NIR_reflection_canopy_coef = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coefficient = Coefficients.Floor.floor_NIR_reflection_coefficient

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

    shScr_PAR_transmission_coefficient = Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coefficient = Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coefficient = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coefficient = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # NIR absorption coefficient of the canopy
    NIR_absorption_canopy_coef = 1 - cover_can_floor_NIR_transmission_coefficient - cover_can_floor_NIR_reflection_coefficient  # page 213
    # NIR absorption coefficient of the floor
    NIR_absorption_floor_coef = cover_can_floor_NIR_transmission_coefficient  # page 213

    cover_PAR_transmission_coefficient = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coefficient, roof_thScr_PAR_transmission_coefficient, shScr_PAR_reflection_coefficient, roof_thScr_PAR_reflection_coefficient)

    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    ratio_GlobPAR = Constants.ratio_GlobPAR
    ratio_GlobNIR = Constants.ratio_GlobNIR
    return ratio_GlobAir * outdoor_global_rad * (cover_PAR_transmission_coefficient * ratio_GlobPAR + (NIR_absorption_canopy_coef + NIR_absorption_floor_coef) * ratio_GlobNIR)


def thermal_screen_FIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.39
    U_ThScr = setpoints.U_ThScr
    thScr_FIR_transmission_coefficient = Coefficients.Thermalscreen.thScr_FIR_transmission_coefficient
    return 1 - U_ThScr * (1 - thScr_FIR_transmission_coefficient)


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


def latent_heat_fluxes(MV) -> float:
    """The latent heat flux
    Equation 8.42
    :param float MV: the vapor flux from object 1 to object 2
    :return: The latent heat flux from object 1 to object 2 [W m^-2]
    """
    return Constants.evaporation_latent_heat * MV


def sensible_heat_flux_between_canopy_and_air(states: States):
    can_t = states.can_t
    air_t = states.air_t
    HEC_CanAir = 2 * Constants.canopy_air_convective_heat_exchange_coef * states.leaf_area_index
    return convective_and_conductive_heat_fluxes(HEC_CanAir, can_t, air_t)


def latent_heat_flux_between_canopy_and_air(states: States, setpoints: Setpoints, weather: Weather):
    mass_vapor_flux_CanAir = canopy_transpiration(states, setpoints, weather)
    return latent_heat_fluxes(mass_vapor_flux_CanAir)


def sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints: Setpoints, states: States):
    air_t = states.air_t
    mechcool_t = states.mechcool_t
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    return convective_and_conductive_heat_fluxes(HEC_MechAir, mechcool_t, air_t)


def sensible_heat_flux_between_heating_pipe_and_greenhouse_air(states: States):
    pipe_t = states.pipe_t
    air_t = states.air_t
    pipe_length = Coefficients.Heating.pipe_length
    phi_external_pipe = Coefficients.Heating.phi_external_pipe
    HEC_PipeAir = 1.99 * math.pi * phi_external_pipe * pipe_length * abs(pipe_t - air_t) ** 0.32
    return convective_and_conductive_heat_fluxes(HEC_PipeAir, pipe_t, air_t)


def sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints: Setpoints):
    # Equation 8.53
    U_Blow = setpoints.U_Blow
    heat_cap_Blow = Coefficients.ActiveClimateControl.heat_cap_Blow
    floor_area = Coefficients.Construction.floor_area
    return U_Blow * heat_cap_Blow / floor_area


def sensible_heat_flux_between_buffer_and_greenhouse_air(states: States):
    # Equation 8.57
    HEC_PasAir = Coefficients.ActiveClimateControl.HEC_PasAir
    soil_3_t = states.soil_j_t[2]  # third layer
    air_t = states.air_t
    return HEC_PasAir * (soil_3_t - air_t)


def sensible_heat_flux_between_floor_and_greenhouse_air(states: States):
    floor_t = states.floor_t
    air_t = states.air_t
    if floor_t > air_t:
        HEC_AirFlr = 1.7 * (floor_t - air_t) ** 0.33
    else:
        HEC_AirFlr = 1.3 * (air_t - floor_t) ** 0.25
    return convective_and_conductive_heat_fluxes(HEC_AirFlr, air_t, floor_t)


def sensible_heat_flux_between_thermal_screen_and_greenhouse_air(states: States, setpoints: Setpoints):
    air_t = states.air_t
    thermal_screen_t = states.thermal_screen_t
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_AirThScr, air_t, thermal_screen_t)


def sensible_heat_flux_between_outdoor_and_greenhouse_air(states: States, setpoints: Setpoints, weather: Weather):
    c_pAir = Constants.c_pAir
    air_t = states.air_t
    outdoor_t = weather.outdoor_t
    density_air = air_density()
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    HEC_AirOut = density_air * c_pAir * (total_side_vent_rate + f_VentForced)
    return convective_and_conductive_heat_fluxes(HEC_AirOut, air_t, outdoor_t)


def sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(states: States, setpoints: Setpoints, weather: Weather):
    c_pAir = Constants.c_pAir
    air_t = states.air_t
    density_air = air_density()
    above_thermal_screen_t = states.above_thermal_screen_t
    thScr_air_flux_rate = thermal_screen_air_flux_rate(setpoints, states, weather)
    HEC_AirTop = density_air * c_pAir * thScr_air_flux_rate
    return convective_and_conductive_heat_fluxes(HEC_AirTop, air_t, above_thermal_screen_t)


def latent_heat_flux_between_fogging_and_greenhouse_air(setpoints: Setpoints):
    mass_vapor_flux_FogAir = fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints)
    return latent_heat_fluxes(mass_vapor_flux_FogAir)


def sensible_heat_flux_between_floor_and_first_layer_soil(states: States):
    floor_thickness = Coefficients.Floor.floor_thickness
    floor_heat_conductivity = Coefficients.Floor.floor_heat_conductivity
    soil_heat_conductivity = Coefficients.Soil.soil_heat_conductivity
    floor_t = states.floor_t
    h_So1 = Coefficients.Soil.soil_thicknesses[0]
    HEC_FlrSo1 = 2 / (floor_thickness / floor_heat_conductivity + h_So1 / soil_heat_conductivity)
    soil_1_t = states.soil_j_t[0]  # first layer
    return convective_and_conductive_heat_fluxes(HEC_FlrSo1, floor_t, soil_1_t)


def latent_heat_flux_between_greenhouse_air_and_thermal_screen(states: States, setpoints: Setpoints):
    air_t = states.air_t
    thermal_screen_t = states.thermal_screen_t
    air_vapor_pressure = saturation_vapor_pressure(air_t)
    thScr_vapor_pressure = saturation_vapor_pressure(thermal_screen_t)
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(air_t - thermal_screen_t) ** 0.33
    mass_vapor_flux_AirThScr = air_to_obj_vapor_flux(air_vapor_pressure, thScr_vapor_pressure, HEC_AirThScr)
    return latent_heat_fluxes(mass_vapor_flux_AirThScr)


def sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(states: States, setpoints: Setpoints):
    above_thermal_screen_t = states.above_thermal_screen_t
    thermal_screen_t = states.thermal_screen_t
    HEC_ThScrTop = 1.7 * setpoints.U_ThScr * abs(thermal_screen_t - above_thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_ThScrTop, thermal_screen_t, above_thermal_screen_t)


def sensible_heat_flux_between_above_thermal_screen_and_internal_cover(states: States, setpoints: Setpoints):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Construction.c_HECin
    cover_area = Coefficients.Construction.cover_area
    floor_area = Coefficients.Construction.floor_area
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * cover_area / floor_area
    return convective_and_conductive_heat_fluxes(HEC_TopCov_in, above_thermal_screen_t, internal_cov_t)


def sensible_heat_flux_between_above_thermal_screen_and_outdoor(states: States, setpoints: Setpoints, weather: Weather):
    c_pAir = Constants.c_pAir
    above_thermal_screen_t = states.above_thermal_screen_t
    outdoor_t = weather.outdoor_t
    density_air = air_density()
    total_roof_vent_rate = total_roof_ventilation_rates(setpoints, states, weather)
    HEC_TopOut = density_air * c_pAir * total_roof_vent_rate
    return convective_and_conductive_heat_fluxes(HEC_TopOut, above_thermal_screen_t, outdoor_t)


def latent_heat_flux_between_above_thermal_screen_and_internal_cover(states: States):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Construction.c_HECin
    cover_area = Coefficients.Construction.cover_area
    floor_area = Coefficients.Construction.floor_area
    above_thermal_screen_vapor_pressure = saturation_vapor_pressure(above_thermal_screen_t)
    cov_in_vapor_pressure = saturation_vapor_pressure(internal_cov_t)
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * cover_area / floor_area
    mass_vapor_flux_TopCov_in = air_to_obj_vapor_flux(above_thermal_screen_vapor_pressure, cov_in_vapor_pressure, HEC_TopCov_in)
    return latent_heat_fluxes(mass_vapor_flux_TopCov_in)


def sensible_heat_flux_between_internal_cover_and_external_cover(states: States, setpoints: Setpoints):
    internal_cov_t = states.internal_cov_t
    external_cov_t = states.external_cov_t
    HEC_Cov_in_Cov_e = lumped_cover_conductive_heat_flux(setpoints)  # Note: line 819 / setGlAux / GreenLight
    return convective_and_conductive_heat_fluxes(HEC_Cov_in_Cov_e, internal_cov_t, external_cov_t)


def sensible_heat_flux_between_external_cover_and_outdoor(states: States, weather: Weather):
    cover_area = Coefficients.Construction.cover_area
    floor_area = Coefficients.Construction.floor_area
    c_HECout_1 = Coefficients.Construction.c_HECout_1
    c_HECout_2 = Coefficients.Construction.c_HECout_2
    external_cov_t = states.external_cov_t
    outdoor_t = weather.outdoor_t
    v_Wind = weather.v_Wind
    c_HECout_3 = Coefficients.Construction.c_HECout_3
    HEC_Cov_e_Out = cover_area * (c_HECout_1 + c_HECout_2 * v_Wind ** c_HECout_3) / floor_area
    return convective_and_conductive_heat_fluxes(HEC_Cov_e_Out, external_cov_t, outdoor_t)


def heat_flux_to_heating_pipe(U, P, A):
    # Equation 8.56
    return U * P / A
