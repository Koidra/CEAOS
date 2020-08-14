"""The implementation of Heat fluxes equations

Based on section 8.6
"""

# 8.6.1 Global, PAR and NIR heat fluxes
from .canopy_transpiration import *
from .radiation_fluxes import *
from .utils import *
from .vapor_fluxes import fogging_system_to_greenhouse_air_latent_vapor_flux, differentiable_air_to_obj_vapor_flux
from ..constants import *


# TODO: inline these
def lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coef):
    # Equation 8.30
    return 1 - cover_NIR_reflection_coef


def floor_virtual_NIR_transmission_coefficients():
    # Equation 8.30
    floor_NIR_reflection_coef = Coefficients.Floor.floor_NIR_reflection_coefficient
    return 1 - floor_NIR_reflection_coef


def canopy_virtual_NIR_transmission_coefficient(states: ClimateStates):
    # Equation 8.31
    return math.exp(-CANOPY_NIR_EXTINCTION_COEF * states.leaf_area_index)


def canopy_virtual_NIR_reflection_coefficient(states: ClimateStates):
    # Equation 8.32
    virtual_NIR_transmission_canopy_coef = canopy_virtual_NIR_transmission_coefficient(states)
    return CANOPY_NIR_REFLECTION_COEF * (1 - virtual_NIR_transmission_canopy_coef)


def greenhouse_air_absorbed_global_radiation(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    """
    Equation 8.36
    Args:
        states:
        setpoints:
        weather:

    Returns: Global radiation from the sun absorbed by the greenhouse air [W m^-2]

    """
    cover_canopy_floor_NIR_transmission_coef, cover_canopy_floor_NIR_reflection_coef = cover_canopy_floor_NIR_coefficients(states, setpoints)

    # NIR absorption coefficient of the canopy
    NIR_absorption_canopy_coef = absorption_coefficient(cover_canopy_floor_NIR_transmission_coef, cover_canopy_floor_NIR_reflection_coef)  # page 213
    # NIR absorption coefficient of the floor
    NIR_absorption_floor_coef = cover_canopy_floor_NIR_transmission_coef  # page 213

    cover_PAR_transmission_coef, cover_PAR_reflection_coef = lumped_cover_PAR_coefficients(setpoints)

    return Coefficients.Construction.ratio_GlobAir * weather.outdoor_global_rad \
           * (cover_PAR_transmission_coef * RATIO_GLOBALPAR
             + (NIR_absorption_canopy_coef + NIR_absorption_floor_coef) * RATIO_GLOBALNIR)


def lamp_radiation(states, setpoints):
    """
    Equation A23 [5]
    Args:
        states:
        setpoints:

    Returns: PAR and NIR from the lamps absorbed by the greenhouse air [W m^{-2}]

    """
    electrical_input_lamp = lamp_electrical_input(setpoints)
    radiation_flux_PAR_LampCanopy = canopy_PAR_absorbed_from_lamp(states, setpoints)
    radiation_flux_NIR_LampCanopy = canopy_NIR_absorbed_from_lamp(states, setpoints)
    radiation_flux_PAR_LampFlr = floor_PAR_absorbed_from_lamp(states, setpoints)
    radiation_flux_NIR_LampFlr = floor_NIR_absorbed_from_lamp(states, setpoints)

    return (Coefficients.Lamp.lamp_electrical_input_PAR_conversion + Coefficients.Lamp.lamp_electrical_input_NIR_conversion) \
           * electrical_input_lamp \
           - radiation_flux_PAR_LampCanopy - radiation_flux_NIR_LampCanopy - radiation_flux_PAR_LampFlr - radiation_flux_NIR_LampFlr


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


def sensible_heat_flux_between_canopy_and_air(states: ClimateStates):
    HEC_CanopyAir = 2 * CANOPY_AIR_CONVECTIVE_HEAT_EXCHANGE_COEF * states.leaf_area_index
    return convective_and_conductive_heat_fluxes(HEC_CanopyAir, states.t_Canopy, states.t_Air)


def latent_heat_flux_between_canopy_and_air(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    mass_vapor_flux_CanopyAir = canopy_transpiration(states, setpoints, weather)
    return latent_heat_fluxes(mass_vapor_flux_CanopyAir)


def sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints: Setpoints, states: ClimateStates):
    HEC_MechAir = mechanical_cooling_to_greenhouse_air_heat_exchange_coefficient(setpoints, states)
    return convective_and_conductive_heat_fluxes(HEC_MechAir, states.t_MechCool, states.t_Air)


def sensible_heat_flux_between_heating_pipe_and_greenhouse_air(states: ClimateStates):
    HEC_PipeAir = 1.99 * math.pi * Coefficients.Heating.phi_external_pipe * Coefficients.Heating.pipe_length \
                  * abs(states.t_Pipe - states.t_Air) ** 0.32
    return convective_and_conductive_heat_fluxes(HEC_PipeAir, states.t_Pipe, states.t_Air)


def sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints: Setpoints):
    # Equation 8.53
    return setpoints.U_Blow * Coefficients.ActiveClimateControl.heat_cap_Blow / Coefficients.Construction.floor_area


def sensible_heat_flux_between_buffer_and_greenhouse_air(states: ClimateStates):
    # Equation 8.57
    soil_3_t = states.t_Soil[2]  # third layer
    return Coefficients.ActiveClimateControl.HEC_PasAir * (soil_3_t - states.t_Air)


def sensible_heat_flux_between_floor_and_greenhouse_air(states: ClimateStates):
    if states.t_Floor > states.t_Air:
        HEC_AirFlr = 1.7 * (states.t_Floor - states.t_Air) ** 0.33
    else:
        HEC_AirFlr = 1.3 * (states.t_Air - states.t_Floor) ** 0.25
    return convective_and_conductive_heat_fluxes(HEC_AirFlr, states.t_Air, states.t_Floor)


def sensible_heat_flux_between_thermal_screen_and_greenhouse_air(states: ClimateStates, setpoints: Setpoints):
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(states.t_Air - states.t_ThScr) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_AirThScr, states.t_Air, states.t_ThScr)


def sensible_heat_flux_between_outdoor_and_greenhouse_air(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    density_air = air_density()
    total_side_vent_rate = total_side_vents_ventilation_rates(setpoints, states, weather)
    f_VentForced = 0  # According to GreenLight, forced ventilation doesn't exist in this greenhouse
    HEC_AirOut = density_air * C_PAIR * (total_side_vent_rate + f_VentForced)
    return convective_and_conductive_heat_fluxes(HEC_AirOut, states.t_Air, weather.t_Outdoor)


def sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    density_air = air_density()
    thScr_air_flux_rate = thermal_screen_air_flux_rate(setpoints, states, weather)
    HEC_AirTop = density_air * C_PAIR * thScr_air_flux_rate
    return convective_and_conductive_heat_fluxes(HEC_AirTop, states.t_Air, states.t_AboveThScr)


def latent_heat_flux_between_fogging_and_greenhouse_air(setpoints: Setpoints):
    mass_vapor_flux_FogAir = fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints)
    return latent_heat_fluxes(mass_vapor_flux_FogAir)


def sensible_heat_flux_between_floor_and_first_layer_soil(states: ClimateStates):
    floor_thickness = Coefficients.Floor.floor_thickness
    floor_heat_conductivity = Coefficients.Floor.floor_heat_conductivity
    soil_heat_conductivity = Coefficients.Soil.soil_heat_conductivity
    h_So1 = Coefficients.Soil.soil_thicknesses[0]
    HEC_FlrSo1 = 2 / (floor_thickness / floor_heat_conductivity + h_So1 / soil_heat_conductivity)
    soil_1_t = states.t_Soil[0]  # first layer
    return convective_and_conductive_heat_fluxes(HEC_FlrSo1, states.t_Floor, soil_1_t)


def latent_heat_flux_between_greenhouse_air_and_thermal_screen(states: ClimateStates, setpoints: Setpoints):
    vapor_pressure_ThScr = saturation_vapor_pressure(states.t_ThScr)
    HEC_AirThScr = 1.7 * setpoints.U_ThScr * abs(states.t_Air - states.t_ThScr) ** 0.33
    mass_vapor_flux_AirThScr = differentiable_air_to_obj_vapor_flux(states.vapor_pressure_Air,
                                                                    vapor_pressure_ThScr,
                                                                    HEC_AirThScr)
    return latent_heat_fluxes(mass_vapor_flux_AirThScr)


def sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(states: ClimateStates, setpoints: Setpoints):
    HEC_ThScrTop = 1.7 * setpoints.U_ThScr * abs(states.t_ThScr - states.t_AboveThScr) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_ThScrTop, states.t_ThScr, states.t_AboveThScr)


def sensible_heat_flux_between_above_thermal_screen_and_internal_cover(states: ClimateStates):
    HEC_TopCov_in = Coefficients.Construction.c_HECin * (states.t_AboveThScr - states.t_Cov_internal) ** 0.33 \
                    * Coefficients.Construction.cover_area / Coefficients.Construction.floor_area
    return convective_and_conductive_heat_fluxes(HEC_TopCov_in, states.t_AboveThScr, states.t_Cov_internal)


def sensible_heat_flux_between_above_thermal_screen_and_outdoor(states: ClimateStates, setpoints: Setpoints, weather: Weather):
    density_air = air_density()
    total_roof_vent_rate = total_roof_ventilation_rates(setpoints, states, weather)
    HEC_TopOut = density_air * C_PAIR * total_roof_vent_rate
    return convective_and_conductive_heat_fluxes(HEC_TopOut, states.t_AboveThScr, weather.t_Outdoor)


def latent_heat_flux_between_above_thermal_screen_and_internal_cover(states: ClimateStates):
    vapor_pressure_Cov_internal = saturation_vapor_pressure(states.t_Cov_internal)
    HEC_TopCov_in = Coefficients.Construction.c_HECin * (states.t_AboveThScr - states.t_Cov_internal) ** 0.33 \
                    * Coefficients.Construction.cover_area / Coefficients.Construction.floor_area
    mass_vapor_flux_TopCov_in = differentiable_air_to_obj_vapor_flux(states.vapor_pressure_AboveThScr,
                                                                     vapor_pressure_Cov_internal,
                                                                     HEC_TopCov_in)
    return latent_heat_fluxes(mass_vapor_flux_TopCov_in)


def sensible_heat_flux_between_internal_cover_and_external_cover(states: ClimateStates):
    HEC_Cov_in_Cov_e = lumped_cover_conductive_heat_flux()  # Note: line 819 / setGlAux / GreenLight
    return convective_and_conductive_heat_fluxes(HEC_Cov_in_Cov_e, states.t_Cov_internal, states.t_Cov_external)


def sensible_heat_flux_between_external_cover_and_outdoor(states: ClimateStates, weather: Weather):
    HEC_Cov_e_Out = Coefficients.Construction.cover_area \
                    * (Coefficients.Construction.c_HECout_1
                       + Coefficients.Construction.c_HECout_2 * weather.v_Wind ** Coefficients.Construction.c_HECout_3) \
                    / Coefficients.Construction.floor_area
    return convective_and_conductive_heat_fluxes(HEC_Cov_e_Out, states.t_Cov_external, weather.t_Outdoor)


def heat_flux_to_heating_pipe(U, P, A):
    # Equation 8.56
    return U * P / A


def sensible_heat_flux_between_greenhouse_air_and_blackout_screen(states: ClimateStates, setpoints: Setpoints):
    """
    Equations A28, A32 [2]
    Args:
        setpoints:
        states:
    Returns: Between air in main compartment and blackout screen [W m^{-2}]
    """
    HEC_AirBlScr = 1.7 * setpoints.U_BlScr * abs(states.t_Air - states.t_BlScr) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_AirBlScr, states.t_Air, states.t_BlScr)


def sensible_heat_flux_between_lamps_and_greenhouse_air(states: ClimateStates):
    """
    Equation A29 [5]
    Args:
        states:

    Returns: Between lamps and air in main compartment [W m^{-2}]
    """
    return convective_and_conductive_heat_fluxes(Coefficients.Lamp.c_HEC_LampAir, states.t_Lamp, states.t_Air)


def sensible_heat_flux_between_inter_lamp_and_greenhouse_air(states: ClimateStates):
    """
    Equation A30 [5]
    Args:
        states:

    Returns: Between interlights and air in main compartment [W m^{-2}]
    """
    return convective_and_conductive_heat_fluxes(Coefficients.Interlight.c_HEC_InterLampAir, states.t_IntLamp, states.t_Air)


def sensible_heat_flux_between_grow_pipe_and_greenhouse_air(states: ClimateStates):
    """
    Equations A31, A33 [5]
    Args:
        states:

    Returns: Between grow pipes and air in main compartment [W m^{-2}]
    """
    HEC_GroPipeAir = 1.99 * math.pi * Coefficients.GrowPipe.phi_external_pipe * Coefficients.GrowPipe.pipe_length * (abs(states.t_GrowPipe - states.t_Air)) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_GroPipeAir, states.t_GrowPipe, states.t_Air)


def sensible_heat_flux_between_above_thermal_screen_and_blackout_screen(states: ClimateStates, setpoints: Setpoints):
    HEC_ThScrBlScr = 1.7 * setpoints.U_BlScr * abs(states.t_BlScr - states.t_AboveThScr) ** 0.33
    return convective_and_conductive_heat_fluxes(HEC_ThScrBlScr, states.t_BlScr, states.t_AboveThScr)


def latent_heat_flux_between_greenhouse_air_and_blackout_screen(states: ClimateStates, setpoints: Setpoints):
    blScr_vapor_pressure = saturation_vapor_pressure(states.t_BlScr)
    HEC_AirBlScr = 1.7 * setpoints.U_BlScr * abs(states.t_Air - states.t_BlScr) ** 0.33
    """
    % Condensation from main compartment on blackout screen [kg m^{-2} s^{-1}]
    % Equation A39 [2]
    """
    mass_vapor_flux_AirBlScr = differentiable_air_to_obj_vapor_flux(states.vapor_pressure_Air,
                                                                    blScr_vapor_pressure, HEC_AirBlScr)
    return latent_heat_fluxes(mass_vapor_flux_AirBlScr)


def sensible_heat_flux_between_boiler_and_grow_pipe(setpoints: Setpoints):
    return setpoints.U_BoilGro*Coefficients.GrowPipe.cap_BoilGro/Coefficients.Construction.floor_area
