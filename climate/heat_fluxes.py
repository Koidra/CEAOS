"""The implementation of Heat fluxes equations

Based on section 8.6
"""

# 8.6.1 Global, PAR and NIR heat fluxes
from climate.canopy_transpiration import *
from climate.electrical_input import lamp_electrical_input
from climate.radiation_fluxes import *
from climate.utils import *
from climate.vapor_fluxes import fogging_system_to_greenhouse_air_latent_vapor_flux, differentiable_air_to_obj_vapor_flux
from constants import *


# TODO: inline these
def lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coef):
    # Equation 8.30
    return 1 - cover_NIR_reflection_coef


def floor_virtual_NIR_transmission_coefficients():
    # Equation 8.30
    floor_NIR_reflection_coef = Coefficients.Floor.floor_NIR_reflection_coefficient
    return 1 - floor_NIR_reflection_coef


def canopy_virtual_NIR_transmission_coefficient(states: States):
    # Equation 8.31
    return math.exp(-CANOPY_NIR_EXTINCTION_COEF * states.leaf_area_index)


def canopy_virtual_NIR_reflection_coefficient(states: States):
    # Equation 8.32
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    return CANOPY_NIR_REFLECTION_COEF * (1 - virtual_NIR_transmission_canopy_coef)


def construction_elements_global_radiation(states: States, setpoints: Setpoints, weather: Weather):
    """
    Equation 8.36
    Args:
        states:
        setpoints:
        weather:

    Returns:[W m^-2]

    """
    outdoor_global_rad = weather.outdoor_global_rad
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

    shScr_PAR_transmission_coef = Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient  # line 156 / setGlParams / GreenLight
    shScr_PAR_reflection_coef = Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 153 / setGlParams / GreenLight

    # PAR transmission coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_transmission_coef = roof_thermal_screen_PAR_transmission_coefficient(setpoints)
    # PAR reflection coefficient of the movable shading screen and the semi-permanent shading screen
    roof_thScr_PAR_reflection_coef = roof_thermal_screen_PAR_reflection_coefficient(setpoints)

    # NIR absorption coefficient of the canopy
    NIR_absorption_canopy_coef = 1 - cover_canopy_floor_NIR_transmission_coef - cover_canopy_floor_NIR_reflection_coef  # page 213
    # NIR absorption coefficient of the floor
    NIR_absorption_floor_coef = cover_canopy_floor_NIR_transmission_coef  # page 213

    cover_PAR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coef,
                                                                              roof_thScr_PAR_transmission_coef,
                                                                              shScr_PAR_reflection_coef,
                                                                              roof_thScr_PAR_reflection_coef)

    ratio_GlobAir = Coefficients.Construction.ratio_GlobAir
    return ratio_GlobAir * outdoor_global_rad * \
           (cover_PAR_transmission_coef * RATIO_GLOBALPAR + (NIR_absorption_canopy_coef + NIR_absorption_floor_coef) * RATIO_GLOBALNIR)


def lamp_radiation(states, setpoints, weather):
    """
    Equation A23 [5]
    Args:
        states:
        setpoints:
        weather:

    Returns: PAR and NIR from the lamps absorbed by the greenhouse air [W m^{-2}]

    """
    electrical_input_lamp = lamp_electrical_input(setpoints)
    radiation_flux_PAR_LampCanopy = canopy_PAR_absorbed_from_lamp(states, setpoints)
    radiation_flux_NIR_LampCanopy = canopy_NIR_absorbed_from_lamp(states, setpoints)
    radiation_flux_PAR_LampFlr = floor_PAR_absorbed_from_lamp(states, setpoints, weather)
    radiation_flux_NIR_LampFlr = floor_NIR_absorbed_from_lamp(states, setpoints)

    return (Coefficients.Lamp.lamp_electrical_input_PAR_conversion + Coefficients.Lamp.lamp_electrical_input_NIR_conversion) * electrical_input_lamp \
           - radiation_flux_PAR_LampCanopy - radiation_flux_NIR_LampCanopy - radiation_flux_PAR_LampFlr - radiation_flux_NIR_LampFlr


def thermal_screen_FIR_transmission_coefficient(setpoints: Setpoints):
    # Equation 8.39
    return 1 - setpoints.U_ThScr * (1 - Coefficients.Thermalscreen.thScr_FIR_transmission_coefficient)


def blackout_screen_FIR_transmission_coefficient(setpoints: Setpoints):
    return 1 - setpoints.U_BlScr * (1 - Coefficients.Blackoutscreen.blScr_FIR_transmission_coef)


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
    return EVAPORATION_LATENT_HEAT * MV


def sensible_heat_flux_between_canopy_and_air(states: States):
    HEC_CanopyAir = 2 * CANOPY_AIR_CONVECTIVE_HEAT_EXCHANGE_COEF * states.leaf_area_index
    return convective_and_conductive_heat_fluxes(HEC_CanopyAir, states.canopy_t, states.air_t)


def latent_heat_flux_between_canopy_and_air(states: States, setpoints: Setpoints, weather: Weather):
    mass_vapor_flux_CanopyAir = canopy_transpiration(states, setpoints, weather)
    return latent_heat_fluxes(mass_vapor_flux_CanopyAir)


def sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints: Setpoints, states: States):
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    return convective_and_conductive_heat_fluxes(HEC_MechAir, states.mechcool_t, states.air_t)


def sensible_heat_flux_between_heating_pipe_and_greenhouse_air(states: States):
    HEC_PipeAir = 1.99 * math.pi * Coefficients.Heating.phi_external_pipe * Coefficients.Heating.pipe_length \
                  * abs(states.pipe_t - states.air_t) ** 0.32
    return convective_and_conductive_heat_fluxes(HEC_PipeAir, states.pipe_t, states.air_t)


def sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints: Setpoints):
    # Equation 8.53
    return setpoints.U_Blow * Coefficients.ActiveClimateControl.heat_cap_Blow / Coefficients.Construction.floor_area


def sensible_heat_flux_between_buffer_and_greenhouse_air(states: States):
    # Equation 8.57
    soil_3_t = states.soil_j_t[2]  # third layer
    return Coefficients.ActiveClimateControl.HEC_PasAir * (soil_3_t - states.air_t)


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
    density_air = air_density()
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    HEC_AirOut = density_air * C_PAIR * (total_side_vent_rate + f_VentForced)
    return convective_and_conductive_heat_fluxes(HEC_AirOut, states.air_t, weather.outdoor_t)


def sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(states: States, setpoints: Setpoints, weather: Weather):
    density_air = air_density()
    thScr_air_flux_rate = thermal_screen_air_flux_rate(setpoints, states, weather)
    HEC_AirTop = density_air * C_PAIR * thScr_air_flux_rate
    return convective_and_conductive_heat_fluxes(HEC_AirTop, states.air_t, states.above_thermal_screen_t)


def latent_heat_flux_between_fogging_and_greenhouse_air(setpoints: Setpoints):
    mass_vapor_flux_FogAir = fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints)
    return latent_heat_fluxes(mass_vapor_flux_FogAir)


def sensible_heat_flux_between_floor_and_first_layer_soil(states: States):
    floor_thickness = Coefficients.Floor.floor_thickness
    floor_heat_conductivity = Coefficients.Floor.floor_heat_conductivity
    soil_heat_conductivity = Coefficients.Soil.soil_heat_conductivity
    h_So1 = Coefficients.Soil.soil_thicknesses[0]
    HEC_FlrSo1 = 2 / (floor_thickness / floor_heat_conductivity + h_So1 / soil_heat_conductivity)
    soil_1_t = states.soil_j_t[0]  # first layer
    return convective_and_conductive_heat_fluxes(HEC_FlrSo1, states.floor_t, soil_1_t)


def latent_heat_flux_between_greenhouse_air_and_thermal_screen(states: States, setpoints: Setpoints):
    thermal_screen_t = states.thermal_screen_t
    thScr_vp = saturation_vapor_pressure(thermal_screen_t)
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(states.air_t - thermal_screen_t) ** 0.33
    mass_vapor_flux_AirThScr = differentiable_air_to_obj_vapor_flux(states.air_vapor_pressure, thScr_vp, HEC_AirThScr)
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
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * Coefficients.Construction.cover_area \
                    / Coefficients.Construction.floor_area
    return convective_and_conductive_heat_fluxes(HEC_TopCov_in, above_thermal_screen_t, internal_cov_t)


def sensible_heat_flux_between_above_thermal_screen_and_outdoor(states: States, setpoints: Setpoints, weather: Weather):
    density_air = air_density()
    total_roof_vent_rate = total_roof_ventilation_rates(setpoints, states, weather)
    HEC_TopOut = density_air * C_PAIR * total_roof_vent_rate
    return convective_and_conductive_heat_fluxes(HEC_TopOut, states.above_thermal_screen_t, weather.outdoor_t)


def latent_heat_flux_between_above_thermal_screen_and_internal_cover(states: States):
    above_thermal_screen_t = states.above_thermal_screen_t
    internal_cov_t = states.internal_cov_t
    c_HECin = Coefficients.Construction.c_HECin
    cov_in_vp = saturation_vapor_pressure(internal_cov_t)
    HEC_TopCov_in = c_HECin * (above_thermal_screen_t - internal_cov_t) ** 0.33 * Coefficients.Construction.cover_area \
                    / Coefficients.Construction.floor_area
    mass_vapor_flux_TopCov_in = differentiable_air_to_obj_vapor_flux(states.above_thermal_screen_vapor_pressure, cov_in_vp, HEC_TopCov_in)
    return latent_heat_fluxes(mass_vapor_flux_TopCov_in)


def sensible_heat_flux_between_internal_cover_and_external_cover(states: States, setpoints: Setpoints):
    HEC_Cov_in_Cov_e = lumped_cover_conductive_heat_flux(setpoints)  # Note: line 819 / setGlAux / GreenLight
    return convective_and_conductive_heat_fluxes(HEC_Cov_in_Cov_e, states.internal_cov_t, states.external_cov_t)


def sensible_heat_flux_between_external_cover_and_outdoor(states: States, weather: Weather):
    HEC_Cov_e_Out = Coefficients.Construction.cover_area \
                    * (Coefficients.Construction.c_HECout_1
                       + Coefficients.Construction.c_HECout_2 * weather.v_Wind ** Coefficients.Construction.c_HECout_3) \
                    / Coefficients.Construction.floor_area
    return convective_and_conductive_heat_fluxes(HEC_Cov_e_Out, states.external_cov_t, weather.outdoor_t)


def heat_flux_to_heating_pipe(U, P, A):
    # Equation 8.56
    return U * P / A


def sensible_heat_flux_between_greenhouse_air_and_blackout_screen(states: States, setpoints: Setpoints):
    """
    Equations A28, A32 [2]
    Args:
        setpoints:
        states:
    Returns: Between air in main compartment and blackout screen [W m^{-2}]
    """
    HEC_AirBlScr = 1.7 * setpoints.U_BlScr * abs(states.air_t - states.blScr_t)**0.33
    return convective_and_conductive_heat_fluxes(HEC_AirBlScr, states.air_t, states.blScr_t)


def sensible_heat_flux_between_lamps_and_greenhouse_air(states: States):
    """
    Equation A29 [5]
    Args:
        states:

    Returns: Between lamps and air in main compartment [W m^{-2}]
    """
    return convective_and_conductive_heat_fluxes(Coefficients.Lamp.c_HEC_LampAir, states.lamp_t, states.air_t)


def sensible_heat_flux_between_inter_lamp_and_greenhouse_air(states: States):
    """
    Equation A30 [5]
    Args:
        states:

    Returns: Between interlights and air in main compartment [W m^{-2}]
    """
    return convective_and_conductive_heat_fluxes(Coefficients.Interlight.c_HEC_InterLampAir, states.intLamp_t, states.air_t)


def sensible_heat_flux_between_grow_pipe_and_greenhouse_air(states: States):
    """
    Equations A31, A33 [5]
    Args:
        states:

    Returns: Between grow pipes and air in main compartment [W m^{-2}]
    """
    HEC_GroPipeAir = 1.99*math.pi*Coefficients.GrowPipe.phi_external_pipe*Coefficients.GrowPipe.pipe_length*(abs(states.groPipe_t-states.air_t))**0.33
    return convective_and_conductive_heat_fluxes(HEC_GroPipeAir, states.groPipe_t, states.air_t)


def sensible_heat_flux_between_above_thermal_screen_and_blackout_screen(states: States, setpoints: Setpoints):
    HEC_ThScrBlScr = 1.7*setpoints.U_BlScr * abs(states.blScr_t - states.above_thermal_screen_t) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_ThScrBlScr, states.blScr_t, states.above_thermal_screen_t)


def latent_heat_flux_between_greenhouse_air_and_blackout_screen(states: States, setpoints: Setpoints):
    blScr_vapor_pressure = saturation_vapor_pressure(states.blScr_t)
    HEC_AirBlScr = 1.7*setpoints.U_BlScr*abs(states.air_t-states.blScr_t)**0.33
    """
    % Condensation from main compartment on blackout screen [kg m^{-2} s^{-1}]
    % Equation A39 [2]
    """
    mass_vapor_flux_AirBlScr = differentiable_air_to_obj_vapor_flux(HEC_AirBlScr,
                                                                    states.air_vapor_pressure,
                                                                    blScr_vapor_pressure)
    return latent_heat_fluxes(mass_vapor_flux_AirBlScr)


def sensible_heat_flux_between_boiler_and_grow_pipe(setpoints: Setpoints):
    return setpoints.U_BoilGro*Coefficients.GrowPipe.cap_BoilGro/Coefficients.Construction.floor_area
