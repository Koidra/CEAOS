"""
The state variables of the model are all described by differential equations.
- can_t: Canopy temperature
- air_t: Greenhouse air temperature
- floor_t: Floor temperature
- soil_j_t: Soil temperature of layer j
- thermal_screen_t: Thermal screen temperature
- above_thermal_screen_t: The air temperature of the compartment above the thermal screen
- internal_cov_t: Internal cover temperature
- external_cov_t: External cover temperature
- pipe_t: Heating pipe temperature
- air_vapor_pressure: Greenhouse air vapor pressure
- top_vapor_pressure: The vapor pressure of the compartment above the thermal screen
- air_CO2: Greenhouse air CO2
- top_CO2: The CO2 of the compartment above the thermal screen
"""
from equations.CO2_fluxes import *
from equations.capacities import *
from equations.heat_fluxes import *
from equations.lumped_cover_layers import lumped_cover_heat_capacity
from equations.radiation_fluxes import *
from equations.vapour_fluxes import *

def canopy_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.1 / 8.1
    # cap_Can * can_t = radiation_flux_PAR_SunCan + radiation_flux_NIR_SunCan + radiation_flux_PipeCan
    #                 - radiation_flux_CanCov_in - radiation_flux_CanFlr - radiation_flux_CanSky - radiation_flux_CanThScr
    #                 - sensible_heat_flux_CanAir - latent_heat_flux_CanAir
    cap_Can = canopy_heat_capacity(inputs)
    radiation_flux_PAR_SunCan = canopy_PAR_absorbed(inputs, setpoints)
    radiation_flux_NIR_SunCan = canopy_NIR_absorbed(inputs, setpoints)
    radiation_flux_PipeCan = FIR_from_pipe_to_canopy(inputs)
    radiation_flux_CanCov_in = FIR_from_canopy_to_internal_cover(inputs, setpoints)
    radiation_flux_CanFlr = FIR_from_canopy_to_floor(inputs)
    radiation_flux_CanSky = FIR_from_canopy_to_sky(inputs, setpoints)
    radiation_flux_CanThScr = FIR_from_canopy_to_thermal_screen(inputs, setpoints)
    sensible_heat_flux_CanAir = sensible_heat_flux_between_canopy_and_air(inputs)
    latent_heat_flux_CanAir = latent_heat_flux_between_canopy_and_air(inputs)

    return (radiation_flux_PAR_SunCan + radiation_flux_NIR_SunCan + radiation_flux_PipeCan
            - radiation_flux_CanCov_in - radiation_flux_CanFlr - radiation_flux_CanSky - radiation_flux_CanThScr
            - sensible_heat_flux_CanAir - latent_heat_flux_CanAir) / cap_Can


def greenhouse_air_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.2 / 8.2
    # cap_Air * air_t = sensible_heat_flux_CanAir + sensible_heat_flux_PadAir + sensible_heat_flux_MechAir
    #                 + sensible_heat_flux_PipeAir + sensible_heat_flux_PasAir + sensible_heat_flux_BlowAir
    #                 + radiation_flux_Glob_SunAir - sensible_heat_flux_AirFlr - sensible_heat_flux_AirThScr
    #                 - sensible_heat_flux_AirOut - sensible_heat_flux_AirTop - sensible_heat_flux_AirOut_Pad
    #                 - latent_heat_flux_AirFog
    h_Air = Constants.Greenhouse.Construction.h_Air
    rho_Air = air_density()
    c_pAir = Constants.Global.c_pAir
    cap_Air = remaining_object_heat_capacity(h_Air, rho_Air, c_pAir)

    radiation_flux_Glob_SunAir = construction_elements_global_radiation(inputs, setpoints)
    sensible_heat_flux_CanAir = sensible_heat_flux_between_canopy_and_air(inputs)
    sensible_heat_flux_PadAir = 0 # sensible_heat_flux_between_pad_and_greenhouse_air(setpoints, inputs)
    sensible_heat_flux_MechAir = sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints, inputs)
    sensible_heat_flux_PipeAir = sensible_heat_flux_between_heating_pipe_and_greenhouse_air(inputs)
    sensible_heat_flux_PasAir = sensible_heat_flux_between_buffer_and_greenhouse_air(inputs)
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    sensible_heat_flux_AirFlr = sensible_heat_flux_between_floor_and_greenhouse_air(inputs)
    sensible_heat_flux_AirThScr = sensible_heat_flux_between_thermal_screen_and_greenhouse_air(inputs, setpoints)
    sensible_heat_flux_AirOut = sensible_heat_flux_between_outdoor_and_greenhouse_air(inputs, setpoints)
    sensible_heat_flux_AirTop = sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(inputs)
    sensible_heat_flux_AirOut_Pad = sensible_heat_flux_between_greenhouse_air_and_outdoor_by_pad_fan_system(setpoints, inputs)
    latent_heat_flux_AirFog = latent_heat_flux_between_fogging_and_greenhouse_air(setpoints)

    return (sensible_heat_flux_CanAir + sensible_heat_flux_PadAir + sensible_heat_flux_MechAir
            + sensible_heat_flux_PipeAir + sensible_heat_flux_PasAir + sensible_heat_flux_BlowAir
            + radiation_flux_Glob_SunAir - sensible_heat_flux_AirFlr - sensible_heat_flux_AirThScr
            - sensible_heat_flux_AirOut - sensible_heat_flux_AirTop - sensible_heat_flux_AirOut_Pad
            - latent_heat_flux_AirFog) / cap_Air


def floor_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.3 / 8.3
    # cap_Flr * floor_t = sensible_heat_flux_AirFlr + radiation_flux_PAR_SunFlr + radiation_flux_NIR_SunFlr
    #                   + radiation_flux_CanFlr + radiation_flux_PipeFlr - sensible_heat_flux_FlrSo1
    #                   - radiation_flux_FlrCov_in - radiation_flux_FlrSky - radiation_flux_FlrThScr
    h_Flr = Constants.Greenhouse.Floor.h_Flr
    rho_Flr = Constants.Greenhouse.Floor.rho_Flr
    c_pFlr = Constants.Greenhouse.Floor.c_pFlr
    cap_Flr = remaining_object_heat_capacity(h_Flr, rho_Flr, c_pFlr)

    sensible_heat_flux_AirFlr = sensible_heat_flux_between_floor_and_greenhouse_air(inputs)
    radiation_flux_PAR_SunFlr = floor_PAR_absorbed(inputs, setpoints)
    radiation_flux_NIR_SunFlr = floor_NIR_absorbed(inputs, setpoints)
    radiation_flux_CanFlr = FIR_from_canopy_to_floor(inputs)
    radiation_flux_PipeFlr = FIR_from_heating_pipe_to_floor(inputs)
    radiation_flux_FlrCov_in = FIR_from_floor_to_internal_cover(inputs, setpoints)
    radiation_flux_FlrSky = FIR_from_floor_to_sky(inputs, setpoints)
    radiation_flux_FlrThScr = FIR_from_floor_to_thermal_screen(inputs, setpoints)
    sensible_heat_flux_FlrSo1 = sensible_heat_flux_between_floor_and_first_layer_soil(inputs)

    return (sensible_heat_flux_AirFlr + radiation_flux_PAR_SunFlr + radiation_flux_NIR_SunFlr
            + radiation_flux_CanFlr + radiation_flux_PipeFlr - sensible_heat_flux_FlrSo1
            - radiation_flux_FlrCov_in - radiation_flux_FlrSky - radiation_flux_FlrThScr) / cap_Flr


def soil_temperature(jth: int, inputs: Inputs): # j = 1,2,..,5
    # Equation 2.4 / 8.4
    # cap_soil_j * soil_j_t = sensible_heat_flux_soil_j_1_soil_j - sensible_heat_flux_soil_j_soil_j_1
    # 0 is Floor, 6 is SoOut

    h_soil_j_minus = Constants.Greenhouse.Floor.h_Flr if jth == 1 else Constants.Greenhouse.Soil.h_So[jth-2]
    h_soil_j = Constants.Greenhouse.Soil.h_So[jth-1]
    h_soil_j_plus = 1.28 if jth == 5 else Constants.Greenhouse.Soil.h_So[jth]# Assumed by GreenLight's authors, line 83, setGlParams
    cap_soil_j = h_soil_j*Constants.Greenhouse.Soil.rho_c_p_So
    lambda_soil = Constants.Greenhouse.Soil.lambda_soil
    HEC_soil_j_1_soil_j = 2 * lambda_soil / (h_soil_j_minus + h_soil_j)
    HEC_soil_j_soil_j_1 = 2 * lambda_soil / (h_soil_j + h_soil_j_plus)
    T_soil_j_minus = inputs.floor_t if jth == 1 else inputs.soil_j_t[jth-2]
    T_soil_j = inputs.soil_j_t[jth-1]
    T_soil_j_plus = inputs.soil_out_t if jth == 5 else inputs.soil_j_t[jth]

    sensible_heat_flux_soil_j_1_soil_j = convective_and_conductive_heat_fluxes(HEC_soil_j_1_soil_j, T_soil_j_minus, T_soil_j)
    sensible_heat_flux_soil_j_soil_j_1 = convective_and_conductive_heat_fluxes(HEC_soil_j_soil_j_1, T_soil_j, T_soil_j_plus)
    return (sensible_heat_flux_soil_j_1_soil_j - sensible_heat_flux_soil_j_soil_j_1) / cap_soil_j


def thermal_screen_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.5 / 8.5
    # cap_ThScr * thermal_screen_t = sensible_heat_flux_AirThScr + latent_heat_flux_AirThScr + radiation_flux_CanThScr
    #                               + radiation_flux_FlrThScr + radiation_flux_PipeThScr - sensible_heat_flux_ThScrTop
    #                               - radiation_flux_ThScrCov_in - radiation_flux_ThScrSky
    h_ThScr = Constants.Greenhouse.Thermalscreen.h_ThScr
    rho_ThScr = Constants.Greenhouse.Thermalscreen.rho_ThScr
    c_pThScr = Constants.Greenhouse.Thermalscreen.c_pThScr
    cap_ThScr = remaining_object_heat_capacity(h_ThScr, rho_ThScr, c_pThScr)

    sensible_heat_flux_AirThScr = sensible_heat_flux_between_thermal_screen_and_greenhouse_air(inputs, setpoints)
    latent_heat_flux_AirThScr = latent_heat_flux_between_greenhouse_air_and_thermal_screen(inputs, setpoints)
    radiation_flux_CanThScr = FIR_from_canopy_to_thermal_screen(inputs, setpoints)
    radiation_flux_FlrThScr = FIR_from_floor_to_thermal_screen(inputs, setpoints)
    radiation_flux_PipeThScr = FIR_from_heating_pipe_to_thermal_screen(inputs, setpoints)
    sensible_heat_flux_ThScrTop = sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(inputs, setpoints)
    radiation_flux_ThScrCov_in = FIR_from_thermal_screen_to_internal_cover(inputs, setpoints)
    radiation_flux_ThScrSky = FIR_from_thermal_screen_to_sky(inputs, setpoints)
    return (sensible_heat_flux_AirThScr + latent_heat_flux_AirThScr + radiation_flux_CanThScr
            + radiation_flux_FlrThScr + radiation_flux_PipeThScr - sensible_heat_flux_ThScrTop
            - radiation_flux_ThScrCov_in - radiation_flux_ThScrSky) / cap_ThScr


def top_compartment_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.6 / 8.6
    # cap_Top * above_thermal_screen_t = sensible_heat_flux_ThScrTop + sensible_heat_flux_AirTop − sensible_heat_flux_TopCov_in − sensible_heat_flux_TopOut
    # TODO: need to recheck if top compartment params are the same with air params
    h_Top = Constants.Greenhouse.Construction.h_Gh - Constants.Greenhouse.Construction.h_Air
    R = Constants.Global.R
    M_Air = Constants.Global.M_Air
    h_Elevation = Constants.Greenhouse.Construction.h_Elevation
    pressure = 101325 * (1 - 2.5577e-5 * h_Elevation) ** 5.25588
    rho_Top = M_Air*pressure/((inputs.above_thermal_screen_t+273.15)*R)  # Note: line 704 / setGlAux / GreenLight
    c_pTop = Constants.Global.c_pAir
    cap_Top = remaining_object_heat_capacity(h_Top, rho_Top, c_pTop)

    sensible_heat_flux_ThScrTop = sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(inputs, setpoints)
    sensible_heat_flux_AirTop = sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(inputs)
    sensible_heat_flux_TopCov_in = sensible_heat_flux_between_above_thermal_screen_and_internal_cover(inputs, setpoints)
    sensible_heat_flux_TopOut = sensible_heat_flux_between_above_thermal_screen_and_outdoor(inputs, setpoints)

    return (sensible_heat_flux_ThScrTop + sensible_heat_flux_AirTop - sensible_heat_flux_TopCov_in - sensible_heat_flux_TopOut) / cap_Top

def internal_cover_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.7 / 8.7
    # cap_Cov_in * internal_cov_t = sensible_heat_flux_TopCov_in + latent_heat_flux_TopCov_in + radiation_flux_CanCov_in
    #                           + radiation_flux_FlrCov_in + radiation_flux_PipeCov_in + radiation_flux_ThScrCov_in
    #                           - sensible_heat_flux_Cov_in_Cov_e
    cap_Cov = lumped_cover_heat_capacity(setpoints)
    cap_Cov_in = internal_external_canopy_heat_capacity(cap_Cov)

    sensible_heat_flux_TopCov_in = sensible_heat_flux_between_above_thermal_screen_and_internal_cover(inputs, setpoints)
    latent_heat_flux_TopCov_in = latent_heat_flux_between_above_thermal_screen_and_internal_cover(inputs)
    radiation_flux_CanCov_in = FIR_from_canopy_to_internal_cover(inputs, setpoints)
    radiation_flux_FlrCov_in = FIR_from_floor_to_internal_cover(inputs, setpoints)
    radiation_flux_PipeCov_in = FIR_from_heating_pipe_to_internal_cover(inputs, setpoints)
    radiation_flux_ThScrCov_in = FIR_from_thermal_screen_to_internal_cover(inputs, setpoints)
    sensible_heat_flux_Cov_in_Cov_e = sensible_heat_flux_between_internal_cover_and_external_cover(inputs, setpoints)
    return (sensible_heat_flux_TopCov_in + latent_heat_flux_TopCov_in + radiation_flux_CanCov_in
            + radiation_flux_FlrCov_in + radiation_flux_PipeCov_in + radiation_flux_ThScrCov_in
            - sensible_heat_flux_Cov_in_Cov_e) / cap_Cov_in


def external_cover_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.8 / 8.8
    # cap_Cov_e * external_cov_t =  radiation_flux_Glob_SunCov_e + sensible_heat_flux_Cov_in_Cov_e
    #                             - sensible_heat_flux_Cov_e_Out - radiation_flux_Cov_e_Sky
    cap_Cov = lumped_cover_heat_capacity(setpoints)
    cap_Cov_e = internal_external_canopy_heat_capacity(cap_Cov)

    radiation_flux_Glob_SunCov_e = cover_global_radiation(inputs, setpoints)
    sensible_heat_flux_Cov_in_Cov_e = sensible_heat_flux_between_internal_cover_and_external_cover(inputs, setpoints)
    sensible_heat_flux_Cov_e_Out = sensible_heat_flux_between_external_cover_and_outdoor(inputs)
    radiation_flux_Cov_e_Sky = FIR_from_external_cover_to_sky(inputs)

    return (radiation_flux_Glob_SunCov_e + sensible_heat_flux_Cov_in_Cov_e - sensible_heat_flux_Cov_e_Out - radiation_flux_Cov_e_Sky) / cap_Cov_e


def heating_pipe_system_surface_temperature(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.9 / 8.9
    # cap_Pipe * pipe_t =  sensible_heat_flux_BoilPipe + sensible_heat_flux_IndPipe + sensible_heat_flux_GeoPipe
    #                    - radiation_flux_PipeSky - radiation_flux_PipeCov_in - radiation_flux_PipeCan
    #                    - radiation_flux_PipeFlr - radiation_flux_PipeThScr - sensible_heat_flux_PipeAir
    cap_Pipe = heating_pipe_heat_capacity()

    U_Boil = setpoints.U_Boil
    P_Boil = Constants.Greenhouse.ActiveClimateControl.P_Boil
    A_Flr = Constants.Greenhouse.Construction.A_Flr
    sensible_heat_flux_BoilPipe = heat_flux_to_heating_pipe(U_Boil, P_Boil, A_Flr)

    U_Ind = setpoints.U_Ind
    P_Ind = Constants.Greenhouse.ActiveClimateControl.P_Ind
    sensible_heat_flux_IndPipe = heat_flux_to_heating_pipe(U_Ind, P_Ind, A_Flr)

    U_Geo = setpoints.U_Geo
    P_Geo = Constants.Greenhouse.ActiveClimateControl.P_Geo
    sensible_heat_flux_GeoPipe = heat_flux_to_heating_pipe(U_Geo, P_Geo, A_Flr)

    radiation_flux_PipeSky = FIR_from_heating_pipe_to_sky(inputs, setpoints)
    radiation_flux_PipeCov_in = FIR_from_heating_pipe_to_internal_cover(inputs, setpoints)
    radiation_flux_PipeCan = FIR_from_pipe_to_canopy(inputs)
    radiation_flux_PipeFlr = FIR_from_heating_pipe_to_floor(inputs)
    radiation_flux_PipeThScr = FIR_from_heating_pipe_to_thermal_screen(inputs, setpoints)
    sensible_heat_flux_PipeAir = sensible_heat_flux_between_heating_pipe_and_greenhouse_air(inputs)

    return (sensible_heat_flux_BoilPipe + sensible_heat_flux_IndPipe + sensible_heat_flux_GeoPipe
            - radiation_flux_PipeSky - radiation_flux_PipeCov_in - radiation_flux_PipeCan
            - radiation_flux_PipeFlr - radiation_flux_PipeThScr - sensible_heat_flux_PipeAir) / cap_Pipe


def greenhouse_air_vapour_pressure(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.10 / 8.10
    # cap_VP_Air * air_vapor_pressure = mass_vapor_flux_CanAir + mass_vapor_flux_PadAir + mass_vapor_flux_FogAir
    #                                 + mass_vapor_flux_BlowAir − mass_vapor_flux_AirThScr − mass_vapor_flux_AirTop
    #                                 − mass_vapor_flux_AirOut − mass_vapor_flux_AirOut_Pad − mass_vapor_flux_AirMech
    cap_VP_Air = air_compartment_water_vapour_capacity(inputs)
    mass_vapor_flux_CanAir = canopy_transpiration(inputs)
    mass_vapor_flux_PadAir = 0 # pad_and_fan_to_greenhouse_air_vapour_flux(setpoints)
    mass_vapor_flux_FogAir = fogging_system_to_greenhouse_air_latent_vapour_flux(setpoints)
    mass_vapor_flux_BlowAir = heat_blower_to_greenhouse_air_vapour_flux(setpoints)
    mass_vapor_flux_AirThScr = greenhouse_air_to_thermal_screen_vapour_flux(setpoints, inputs)
    mass_vapor_flux_AirTop = greenhouse_air_to_above_thermal_screen_vapour_flux(inputs)
    mass_vapor_flux_AirOut = greenhouse_air_to_outdoor_vapour_flux(inputs, setpoints)
    mass_vapor_flux_AirOut_Pad = greenhouse_air_to_outdoor_vapour_flux_by_pad_fan_system(setpoints, inputs)
    mass_vapor_flux_AirMech = greenhouse_air_to_mechanical_cooling_vapour_flux(inputs, setpoints)

    return (mass_vapor_flux_CanAir + mass_vapor_flux_PadAir + mass_vapor_flux_FogAir
            + mass_vapor_flux_BlowAir - mass_vapor_flux_AirThScr - mass_vapor_flux_AirTop
            - mass_vapor_flux_AirOut - mass_vapor_flux_AirOut_Pad - mass_vapor_flux_AirMech) / cap_VP_Air


def top_compartment_vapour_pressure(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.11 / 8.11
    # cap_VP_Top * top_vapor_pressure = mass_vapor_flux_AirTop − mass_vapor_flux_TopCov_in − mass_vapor_flux_TopOut
    cap_VP_Top = air_compartment_water_vapour_capacity(inputs)
    mass_vapor_flux_AirTop = greenhouse_air_to_above_thermal_screen_vapour_flux(inputs)
    mass_vapor_flux_TopCov_in = above_thermal_screen_to_internal_cover_vapour_flux(inputs)
    mass_vapor_flux_TopOut = above_thermal_screen_to_outdoor_vapour_flux(inputs, setpoints)
    return (mass_vapor_flux_AirTop - mass_vapor_flux_TopCov_in - mass_vapor_flux_TopOut) / cap_VP_Top


def greenhouse_air_CO2(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.12 / 8.12
    # cap_CO2_Air * air_CO2 = mass_CO2_flux_BlowAir + mass_CO2_flux_ExtAir + mass_CO2_flux_PadAir
    #                       - mass_CO2_flux_AirCan - mass_CO2_flux_AirTop - mass_CO2_flux_AirOut
    cap_CO2_Air = Constants.Greenhouse.Construction.h_Air # Note: line 45 / setDepParams / GreenLight
    mass_CO2_flux_BlowAir = heat_blower_to_greenhouse_air_CO2_flux(setpoints)
    mass_CO2_flux_ExtAir = external_CO2_added(setpoints)
    mass_CO2_flux_PadAir = 0 # pad_fan_system_and_greenhouse_air_CO2_flux(inputs)
    mass_CO2_flux_AirCan = inputs.MC_AirCan # TODO: check this
    mass_CO2_flux_AirTop = greenhouse_air_and_above_thermal_screen_CO2_flux(inputs)
    mass_CO2_flux_AirOut = greenhouse_air_and_outdoor_CO2_flux(inputs, setpoints)
    return (mass_CO2_flux_BlowAir + mass_CO2_flux_ExtAir + mass_CO2_flux_PadAir
            - mass_CO2_flux_AirCan - mass_CO2_flux_AirTop - mass_CO2_flux_AirOut) / cap_CO2_Air


def top_compartment_air_CO2(setpoints: Setpoints, inputs: Inputs):
    # Equation 2.13 / 8.13
    # cap_CO2_Top * top_CO2 = mass_CO2_flux_AirTop - mass_CO2_flux_TopOut
    cap_CO2_Top = Constants.Greenhouse.Construction.h_Gh - Constants.Greenhouse.Construction.h_Air # Note: line 46 / setDepParams / GreenLight
    mass_CO2_flux_AirTop = greenhouse_air_and_above_thermal_screen_CO2_flux(inputs)
    mass_CO2_flux_TopOut = above_thermal_screen_and_outdoor_CO2_flux(inputs, setpoints)
    return (mass_CO2_flux_AirTop - mass_CO2_flux_TopOut) / cap_CO2_Top
