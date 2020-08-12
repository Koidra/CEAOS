"""
The state variables of the model are all described by differential equations.
- canopy_t: Canopy temperature
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
from .CO2_fluxes import *
from .canopy_transpiration import canopy_transpiration
from .electrical_input import lamp_electrical_input, inter_lamp_electrical_input
from .heat_fluxes import *
from .capacities import *
from .lumped_cover_layers import *
from .radiation_fluxes import *
from .vapor_fluxes import *
from .utils import air_density


def canopy_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.1 / 8.1 [W m-2]
    cap_canopy * canopy_t = radiation_flux_PAR_SunCanopy + radiation_flux_NIR_SunCanopy + radiation_flux_PipeCanopy
                    - radiation_flux_CanopyCov_in - radiation_flux_CanopyFlr - radiation_flux_CanopySky - radiation_flux_CanopyThScr
                    - sensible_heat_flux_CanopyAir - latent_heat_flux_CanopyAir - radiation_flux_CanopyBlScr
                    + radiation_flux_PAR_LampCanopy + radiation_flux_NIR_LampCanopy + radiation_flux_FIR_LampCanopy
                    + radiation_flux_PAR_IntLampCanopy + radiation_flux_NIR_IntLampCanopy + radiation_flux_FIR_IntLampCanopy
                    + radiation_flux_GroPipeCanopy
    :return: The canopy temperature
    """
    cap_canopy = canopy_heat_capacity(states)
    radiation_flux_PAR_SunCanopy = canopy_PAR_absorbed(states, setpoints, weather)
    radiation_flux_NIR_SunCanopy = canopy_NIR_absorbed(states, setpoints, weather)
    radiation_flux_PipeCanopy = FIR_from_pipe_to_canopy(states)
    radiation_flux_CanopyCov_in = FIR_from_canopy_to_internal_cover(states, setpoints)
    radiation_flux_CanopyFlr = FIR_from_canopy_to_floor(states)
    radiation_flux_CanopySky = FIR_from_canopy_to_sky(states, setpoints, weather)
    radiation_flux_CanopyThScr = FIR_from_canopy_to_thermal_screen(states, setpoints)
    sensible_heat_flux_CanopyAir = sensible_heat_flux_between_canopy_and_air(states)
    latent_heat_flux_CanopyAir = latent_heat_flux_between_canopy_and_air(states, setpoints, weather)

    radiation_flux_CanopyBlScr = FIR_from_canopy_to_blackscreen(states, setpoints)
    radiation_flux_PAR_LampCanopy = canopy_PAR_absorbed_from_lamp(states, setpoints)
    radiation_flux_NIR_LampCanopy = canopy_NIR_absorbed_from_lamp(states, setpoints)
    radiation_flux_FIR_LampCanopy = FIR_from_lamp_to_canopy(states)
    radiation_flux_PAR_IntLampCanopy = canopy_PAR_absorbed_from_inter_lamp(setpoints)
    radiation_flux_NIR_IntLampCanopy = canopy_NIR_absorbed_from_inter_lamp(setpoints)
    radiation_flux_FIR_IntLampCanopy = FIR_from_inter_lamp_to_canopy(states)
    radiation_flux_GroPipeCanopy = FIR_from_grow_pipe_to_canopy(states)
    return (radiation_flux_PAR_SunCanopy + radiation_flux_NIR_SunCanopy + radiation_flux_PipeCanopy
            - radiation_flux_CanopyCov_in - radiation_flux_CanopyFlr - radiation_flux_CanopySky - radiation_flux_CanopyThScr
            - sensible_heat_flux_CanopyAir - latent_heat_flux_CanopyAir - radiation_flux_CanopyBlScr
            + radiation_flux_PAR_LampCanopy + radiation_flux_NIR_LampCanopy + radiation_flux_FIR_LampCanopy
            + radiation_flux_PAR_IntLampCanopy + radiation_flux_NIR_IntLampCanopy + radiation_flux_FIR_IntLampCanopy
            + radiation_flux_GroPipeCanopy) / cap_canopy


def greenhouse_air_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.2 / 8.2 [W m-2]
    cap_Air * air_t = sensible_heat_flux_CanopyAir + sensible_heat_flux_MechAir
                    + sensible_heat_flux_PipeAir + sensible_heat_flux_PasAir + sensible_heat_flux_BlowAir
                    + radiation_flux_Glob_SunAir - sensible_heat_flux_AirFlr - sensible_heat_flux_AirThScr
                    - sensible_heat_flux_AirOut - sensible_heat_flux_AirTop
                    - latent_heat_flux_AirFog - sensible_heat_flux_AirBlScr + sensible_heat_flux_LampAir
                    + radiation_flux_LampAir + sensible_heat_flux_IntLampAir + sensible_heat_flux_GroPipeAir
    :return: The greenhouse air temperature
    """
    air_height = coefs.Construction.air_height
    density_air = air_density()
    cap_Air = remaining_object_heat_capacity(air_height, density_air, C_PAIR)

    radiation_flux_Glob_SunAir = construction_elements_global_radiation(states, setpoints, weather)
    sensible_heat_flux_CanopyAir = sensible_heat_flux_between_canopy_and_air(states)
    sensible_heat_flux_MechAir = sensible_heat_flux_between_mechanical_cooling_and_greenhouse_air(setpoints, states)
    sensible_heat_flux_PipeAir = sensible_heat_flux_between_heating_pipe_and_greenhouse_air(states)
    sensible_heat_flux_PasAir = sensible_heat_flux_between_buffer_and_greenhouse_air(states)
    sensible_heat_flux_BlowAir = sensible_heat_flux_between_direct_air_heater_and_greenhouse_air(setpoints)
    sensible_heat_flux_AirFlr = sensible_heat_flux_between_floor_and_greenhouse_air(states)
    sensible_heat_flux_AirThScr = sensible_heat_flux_between_thermal_screen_and_greenhouse_air(states, setpoints)
    sensible_heat_flux_AirOut = sensible_heat_flux_between_outdoor_and_greenhouse_air(states, setpoints, weather)
    sensible_heat_flux_AirTop = sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(states, setpoints, weather)
    latent_heat_flux_AirFog = latent_heat_flux_between_fogging_and_greenhouse_air(setpoints)

    sensible_heat_flux_AirBlScr = sensible_heat_flux_between_greenhouse_air_and_black_screen(states, setpoints)
    sensible_heat_flux_LampAir = sensible_heat_flux_between_lamps_and_greenhouse_air(states)
    radiation_flux_LampAir = lamp_radiation(states, setpoints, weather)
    sensible_heat_flux_IntLampAir = sensible_heat_flux_between_inter_lamp_and_greenhouse_air(states)
    sensible_heat_flux_GroPipeAir = sensible_heat_flux_between_grow_pipe_and_greenhouse_air(states)

    return (sensible_heat_flux_CanopyAir + sensible_heat_flux_MechAir
            + sensible_heat_flux_PipeAir + sensible_heat_flux_PasAir + sensible_heat_flux_BlowAir
            + radiation_flux_Glob_SunAir - sensible_heat_flux_AirFlr - sensible_heat_flux_AirThScr
            - sensible_heat_flux_AirOut - sensible_heat_flux_AirTop - latent_heat_flux_AirFog
            - sensible_heat_flux_AirBlScr + sensible_heat_flux_LampAir + radiation_flux_LampAir
            + sensible_heat_flux_IntLampAir + sensible_heat_flux_GroPipeAir) / cap_Air


def floor_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.3 / 8.3 [W m-2]
    cap_Flr * floor_t = sensible_heat_flux_AirFlr + radiation_flux_PAR_SunFlr + radiation_flux_NIR_SunFlr
                      + radiation_flux_CanopyFlr + radiation_flux_PipeFlr - sensible_heat_flux_FlrSo1
                      - radiation_flux_FlrCov_in - radiation_flux_FlrSky - radiation_flux_FlrThScr
    :return: The floor temperature
    """
    floor_thickness = coefs.Floor.floor_thickness
    floor_density = coefs.Floor.floor_density
    c_pFlr = coefs.Floor.c_pFlr
    cap_Flr = remaining_object_heat_capacity(floor_thickness, floor_density, c_pFlr)

    sensible_heat_flux_AirFlr = sensible_heat_flux_between_floor_and_greenhouse_air(states)
    radiation_flux_PAR_SunFlr = floor_PAR_absorbed(states, setpoints, weather)
    radiation_flux_NIR_SunFlr = floor_NIR_absorbed(states, setpoints, weather)
    radiation_flux_CanopyFlr = FIR_from_canopy_to_floor(states)
    radiation_flux_PipeFlr = FIR_from_heating_pipe_to_floor(states)
    radiation_flux_FlrCov_in = FIR_from_floor_to_internal_cover(states, setpoints)
    radiation_flux_FlrSky = FIR_from_floor_to_sky(states, setpoints, weather)
    radiation_flux_FlrThScr = FIR_from_floor_to_thermal_screen(states, setpoints)
    sensible_heat_flux_FlrSo1 = sensible_heat_flux_between_floor_and_first_layer_soil(states)

    radiation_flux_FlrBlScr = FIR_from_floor_to_black_screen(states, setpoints, weather)
    radiation_flux_PAR_LampFlr = floor_PAR_absorbed_from_lamp(states, setpoints, weather)
    radiation_flux_NIR_LampFlr = floor_NIR_absorbed_from_lamp(states, setpoints, weather)
    radiation_flux_FIR_LampFlr = FIR_from_lamp_to_floor(states, setpoints, weather)

    return (sensible_heat_flux_AirFlr + radiation_flux_PAR_SunFlr + radiation_flux_NIR_SunFlr
            + radiation_flux_CanopyFlr + radiation_flux_PipeFlr - sensible_heat_flux_FlrSo1
            - radiation_flux_FlrCov_in - radiation_flux_FlrSky - radiation_flux_FlrThScr
            - radiation_flux_FlrBlScr
            + radiation_flux_PAR_LampFlr + radiation_flux_NIR_LampFlr + radiation_flux_FIR_LampFlr) / cap_Flr


def soil_temperature(j: int, states: States, weather: Weather):  # j = 1,2,..,5
    """
    Equation 2.4 / 8.4 [W m-2]
    cap_soil_j * soil_j_t = sensible_heat_flux_soil_j_minus_soil_j - sensible_heat_flux_soil_j_soil_j_plus
    0 is Floor, 6 is SoOut
    :return: The soil temperature
    """

    h_soil_j_minus = coefs.Floor.floor_thickness if j == 1 else coefs.Soil.soil_thicknesses[j - 2]
    h_soil_j = coefs.Soil.soil_thicknesses[j - 1]
    h_soil_j_plus = 1.28 if j == 5 else coefs.Soil.soil_thicknesses[j]  # Assumed by GreenLight's authors, line 83, setGlParams
    cap_soil_j = h_soil_j * coefs.Soil.rho_c_p_So
    soil_heat_conductivity = coefs.Soil.soil_heat_conductivity
    HEC_soil_j_minus_soil_j = 2 * soil_heat_conductivity / (h_soil_j_minus + h_soil_j)
    HEC_soil_j_soil_j_plus = 2 * soil_heat_conductivity / (h_soil_j + h_soil_j_plus)
    soil_j_minus_t = states.floor_t if j == 1 else states.soil_j_t[j - 2]
    soil_j_t = states.soil_j_t[j - 1]
    soil_j_plus_t = weather.soil_out_t if j == 5 else states.soil_j_t[j]

    sensible_heat_flux_soil_j_minus_soil_j = convective_and_conductive_heat_fluxes(HEC_soil_j_minus_soil_j, soil_j_minus_t, soil_j_t)
    sensible_heat_flux_soil_j_soil_j_plus = convective_and_conductive_heat_fluxes(HEC_soil_j_soil_j_plus, soil_j_t, soil_j_plus_t)
    return (sensible_heat_flux_soil_j_minus_soil_j - sensible_heat_flux_soil_j_soil_j_plus) / cap_soil_j


def thermal_screen_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.5 / 8.5 [W m-2]
    cap_ThScr * thermal_screen_t = sensible_heat_flux_AirThScr + latent_heat_flux_AirThScr + radiation_flux_CanopyThScr
                                  + radiation_flux_FlrThScr + radiation_flux_PipeThScr - sensible_heat_flux_ThScrTop
                                  - radiation_flux_ThScrCov_in - radiation_flux_ThScrSky + radiation_flux_BlScrThScr
                                  + radiation_flux_LampThScr
    :return: The thermal screen temperature
    """
    thScr_thickness = coefs.Thermalscreen.thScr_thickness
    thScr_density = coefs.Thermalscreen.thScr_density
    c_pThScr = coefs.Thermalscreen.c_pThScr
    cap_ThScr = remaining_object_heat_capacity(thScr_thickness, thScr_density, c_pThScr)

    sensible_heat_flux_AirThScr = sensible_heat_flux_between_thermal_screen_and_greenhouse_air(states, setpoints)
    latent_heat_flux_AirThScr = latent_heat_flux_between_greenhouse_air_and_thermal_screen(states, setpoints)
    radiation_flux_CanopyThScr = FIR_from_canopy_to_thermal_screen(states, setpoints)
    radiation_flux_FlrThScr = FIR_from_floor_to_thermal_screen(states, setpoints)
    radiation_flux_PipeThScr = FIR_from_heating_pipe_to_thermal_screen(states, setpoints)
    sensible_heat_flux_ThScrTop = sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(states, setpoints)
    radiation_flux_ThScrCov_in = FIR_from_thermal_screen_to_internal_cover(states, setpoints)
    radiation_flux_ThScrSky = FIR_from_thermal_screen_to_sky(states, setpoints, weather)
    radiation_flux_BlScrThScr = FIR_from_black_screen_to_thermal_screen(states, setpoints, weather)
    radiation_flux_LampThScr = FIR_from_lamp_to_thermal_screen(states, setpoints, weather)
    return (sensible_heat_flux_AirThScr + latent_heat_flux_AirThScr + radiation_flux_CanopyThScr
            + radiation_flux_FlrThScr + radiation_flux_PipeThScr - sensible_heat_flux_ThScrTop
            - radiation_flux_ThScrCov_in - radiation_flux_ThScrSky + radiation_flux_BlScrThScr
            + radiation_flux_LampThScr) / cap_ThScr


def top_compartment_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.6 / 8.6 [W m-2]
    cap_Top * above_thermal_screen_t = sensible_heat_flux_ThScrTop + sensible_heat_flux_AirTop − sensible_heat_flux_TopCov_in
                                     − sensible_heat_flux_TopOut + sensible_heat_flux_BlScrTop
    TODO: need to recheck if top compartment params are the same with air params
    :return: The above thermal screen air temperature
    """
    h_Top = coefs.Construction.greenhouse_height - coefs.Construction.air_height
    elevation_height = coefs.Construction.elevation_height
    pressure = 101325 * (1 - 2.5577e-5 * elevation_height) ** 5.25588
    density_Top = M_AIR * pressure / ((states.above_thermal_screen_t + 273.15) * M_GAS)  # Note: line 704 / setGlAux / GreenLight
    c_pTop = C_PAIR
    cap_Top = remaining_object_heat_capacity(h_Top, density_Top, c_pTop)

    sensible_heat_flux_ThScrTop = sensible_heat_flux_between_thermal_screen_and_above_thermal_screen(states, setpoints)
    sensible_heat_flux_AirTop = sensible_heat_flux_between_above_thermal_screen_and_greenhouse_air(states, setpoints, weather)
    sensible_heat_flux_TopCov_in = sensible_heat_flux_between_above_thermal_screen_and_internal_cover(states, setpoints)
    sensible_heat_flux_TopOut = sensible_heat_flux_between_above_thermal_screen_and_outdoor(states, setpoints, weather)
    sensible_heat_flux_BlScrTop = sensible_heat_flux_between_above_thermal_screen_and_black_screen(states)
    return (sensible_heat_flux_ThScrTop + sensible_heat_flux_AirTop
            - sensible_heat_flux_TopCov_in - sensible_heat_flux_TopOut + sensible_heat_flux_BlScrTop) / cap_Top


def internal_cover_temperature(setpoints: Setpoints, states: States):
    """
    Equation 2.7 / 8.7 [W m-2]
    cap_Cov_in * internal_cov_t = sensible_heat_flux_TopCov_in + latent_heat_flux_TopCov_in + radiation_flux_CanopyCov_in
                              + radiation_flux_FlrCov_in + radiation_flux_PipeCov_in + radiation_flux_ThScrCov_in
                              - sensible_heat_flux_Cov_in_Cov_e + radiation_flux_BlScrCov_in + radiation_flux_LampCov_in
    :return: The internal cover temperature
    """
    cap_Cov = lumped_cover_heat_capacity(setpoints)
    cap_Cov_in = internal_external_canopy_heat_capacity(cap_Cov)

    sensible_heat_flux_TopCov_in = sensible_heat_flux_between_above_thermal_screen_and_internal_cover(states, setpoints)
    latent_heat_flux_TopCov_in = latent_heat_flux_between_above_thermal_screen_and_internal_cover(states)
    radiation_flux_CanopyCov_in = FIR_from_canopy_to_internal_cover(states, setpoints)
    radiation_flux_FlrCov_in = FIR_from_floor_to_internal_cover(states, setpoints)
    radiation_flux_PipeCov_in = FIR_from_heating_pipe_to_internal_cover(states, setpoints)
    radiation_flux_ThScrCov_in = FIR_from_thermal_screen_to_internal_cover(states, setpoints)
    sensible_heat_flux_Cov_in_Cov_e = sensible_heat_flux_between_internal_cover_and_external_cover(states, setpoints)
    radiation_flux_BlScrCov_in = FIR_from_black_screen_to_internal_cover(states, setpoints)
    radiation_flux_LampCov_in = FIR_from_lamp_to_internal_cover(states, setpoints)
    return (sensible_heat_flux_TopCov_in + latent_heat_flux_TopCov_in + radiation_flux_CanopyCov_in
            + radiation_flux_FlrCov_in + radiation_flux_PipeCov_in + radiation_flux_ThScrCov_in
            - sensible_heat_flux_Cov_in_Cov_e + radiation_flux_BlScrCov_in + radiation_flux_LampCov_in) / cap_Cov_in


def external_cover_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.8 / 8.8 [W m-2]
    cap_Cov_e * external_cov_t =  radiation_flux_Glob_SunCov_e + sensible_heat_flux_Cov_in_Cov_e
                                - sensible_heat_flux_Cov_e_Out - radiation_flux_Cov_e_Sky
    :return: The external cover temperature
    """
    cap_Cov = lumped_cover_heat_capacity(setpoints)
    cap_Cov_e = internal_external_canopy_heat_capacity(cap_Cov)

    radiation_flux_Glob_SunCov_e = cover_global_radiation(states, setpoints, weather)
    sensible_heat_flux_Cov_in_Cov_e = sensible_heat_flux_between_internal_cover_and_external_cover(states, setpoints)
    sensible_heat_flux_Cov_e_Out = sensible_heat_flux_between_external_cover_and_outdoor(states, weather)
    radiation_flux_Cov_e_Sky = FIR_from_external_cover_to_sky(states, weather)

    return (radiation_flux_Glob_SunCov_e + sensible_heat_flux_Cov_in_Cov_e
            - sensible_heat_flux_Cov_e_Out - radiation_flux_Cov_e_Sky) / cap_Cov_e


def heating_pipe_system_surface_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.9 / 8.9 [W m-2]
    cap_Pipe * pipe_t =  sensible_heat_flux_BoilPipe + sensible_heat_flux_IndPipe + sensible_heat_flux_GeoPipe
                       - radiation_flux_PipeSky - radiation_flux_PipeCov_in - radiation_flux_PipeCanopy
                       - radiation_flux_PipeFlr - radiation_flux_PipeThScr - sensible_heat_flux_PipeAir
                       - radiation_flux_PipeBlScr + radiation_flux_LampPipe
    :return: The heating pipe temperature
    """
    cap_Pipe = heating_pipe_heat_capacity()

    U_Boil = setpoints.U_Boil
    heat_cap_Boil = Coefficients.ActiveClimateControl.heat_cap_Boil
    floor_area = Coefficients.Construction.floor_area
    sensible_heat_flux_BoilPipe = heat_flux_to_heating_pipe(U_Boil, heat_cap_Boil, floor_area)

    U_Ind = setpoints.U_Ind
    heat_cap_Ind = Coefficients.ActiveClimateControl.heat_cap_Ind
    sensible_heat_flux_IndPipe = heat_flux_to_heating_pipe(U_Ind, heat_cap_Ind, floor_area)

    U_Geo = setpoints.U_Geo
    heat_cap_Geo = Coefficients.ActiveClimateControl.heat_cap_Geo
    sensible_heat_flux_GeoPipe = heat_flux_to_heating_pipe(U_Geo, heat_cap_Geo, floor_area)

    radiation_flux_PipeSky = FIR_from_heating_pipe_to_sky(states, setpoints, weather)
    radiation_flux_PipeCov_in = FIR_from_heating_pipe_to_internal_cover(states, setpoints)
    radiation_flux_PipeCanopy = FIR_from_pipe_to_canopy(states)
    radiation_flux_PipeFlr = FIR_from_heating_pipe_to_floor(states)
    radiation_flux_PipeThScr = FIR_from_heating_pipe_to_thermal_screen(states, setpoints)
    sensible_heat_flux_PipeAir = sensible_heat_flux_between_heating_pipe_and_greenhouse_air(states)
    radiation_flux_PipeBlScr = FIR_from_heating_pipe_to_black_screen(states, setpoints, weather)
    radiation_flux_LampPipe = FIR_from_lamp_to_heating_pipe(states, setpoints, weather)
    return (sensible_heat_flux_BoilPipe + sensible_heat_flux_IndPipe + sensible_heat_flux_GeoPipe
            - radiation_flux_PipeSky - radiation_flux_PipeCov_in - radiation_flux_PipeCanopy
            - radiation_flux_PipeFlr - radiation_flux_PipeThScr - sensible_heat_flux_PipeAir
            - radiation_flux_PipeBlScr + radiation_flux_LampPipe) / cap_Pipe


def black_screen_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 1 [2] [W m-2]
    cap_BlScr * blScr_t = sensible_heat_flux_AirBlScr + latent_heat_flux_AirBlScr + radiation_flux_CanBlScr
                        + radiation_flux_FlrBlScr + radiation_flux_PipeBlScr - sensible_heat_flux_BlScrTop
                        - radiation_flux_BlScrCov_in - radiation_flux_BlScrSky - radiation_flux_BlScrThScr
                        + radiation_flux_LampBlScr
    Returns: The black screen temperature
    """
    cap_BlScr = 0
    sensible_heat_flux_AirBlScr = sensible_heat_flux_between_greenhouse_air_and_black_screen(states, setpoints)
    latent_heat_flux_AirBlScr = latent_heat_flux_between_greenhouse_air_and_black_screen(states)
    radiation_flux_CanBlScr = sensible_heat_flux_between_canopy_and_black_screen(states)
    radiation_flux_FlrBlScr = FIR_from_floor_to_black_screen(states, setpoints, weather)
    radiation_flux_PipeBlScr = FIR_from_heating_pipe_to_black_screen(states, setpoints, weather)
    sensible_heat_flux_BlScrTop = sensible_heat_flux_between_above_thermal_screen_and_black_screen(states)
    radiation_flux_BlScrCov_in = FIR_from_black_screen_to_internal_cover(states, setpoints)
    radiation_flux_BlScrSky = sensible_heat_flux_between_black_screen_and_sky(states)
    radiation_flux_BlScrThScr = FIR_from_black_screen_to_thermal_screen(states, setpoints, weather)
    radiation_flux_LampBlScr = sensible_heat_flux_between_lamp_and_black_screen(states)

    return (sensible_heat_flux_AirBlScr + latent_heat_flux_AirBlScr + radiation_flux_CanBlScr
            + radiation_flux_FlrBlScr + radiation_flux_PipeBlScr - sensible_heat_flux_BlScrTop
            - radiation_flux_BlScrCov_in - radiation_flux_BlScrSky - radiation_flux_BlScrThScr
            + radiation_flux_LampBlScr)/cap_BlScr


def grow_pipe_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 1 [2] [W m-2]
    cap_GroPipe * groPipe_t = sensible_heat_flux_BoilGroPipe -  radiation_flux_GroPipeCanopy - sensible_heat_flux_GroPipeAir
    Returns: The grow pipe temperature
    """
    cap_GroPipe = 0
    sensible_heat_flux_BoilGroPipe = sensible_heat_flux_between_boiler_and_grow_pipe(states)
    radiation_flux_GroPipeCanopy = FIR_from_grow_pipe_to_canopy(states)
    sensible_heat_flux_GroPipeAir = sensible_heat_flux_between_grow_pipe_and_greenhouse_air(states)
    return (sensible_heat_flux_BoilGroPipe -  radiation_flux_GroPipeCanopy - sensible_heat_flux_GroPipeAir) / cap_GroPipe


def lamps_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2 [2] [W m-2]
    cap_Lamp * lamp_t = electrical_input_Lamp - radiation_flux_LampSky - radiation_flux_LampCov_in
                      - radiation_flux_LampThScr - radiation_flux_LampBlScr - sensible_heat_flux_LampAir
                      - radiation_flux_PAR_LampCanopy - radiation_flux_NIR_LampCanopy - radiation_flux_FIR_LampCanopy
                      - radiation_flux_LampPipe - radiation_flux_PAR_LampFlr - radiation_flux_NIR_LampFlr
                      - radiation_flux_FIR_LampFlr - radiation_flux_LampAir - sensible_heat_flux_LampCool
    Returns:
    """
    cap_Lamp = 0
    electrical_input_Lamp = lamp_electrical_input(setpoints)
    radiation_flux_LampSky = FIR_from_lamp_to_sky(states, setpoints, weather)
    radiation_flux_LampCov_in = FIR_from_lamp_to_internal_cover(states, setpoints)
    radiation_flux_LampThScr = FIR_from_lamp_to_thermal_screen(states, setpoints, weather)
    radiation_flux_LampBlScr = sensible_heat_flux_between_lamp_and_black_screen(states)
    sensible_heat_flux_LampAir = sensible_heat_flux_between_lamps_and_greenhouse_air(states)
    radiation_flux_PAR_LampCanopy = canopy_PAR_absorbed_from_lamp(states, setpoints)
    radiation_flux_NIR_LampCanopy = canopy_NIR_absorbed_from_lamp(states, setpoints)
    radiation_flux_FIR_LampCanopy = FIR_from_inter_lamp_to_canopy(states)
    radiation_flux_LampPipe = FIR_from_lamp_to_heating_pipe(states, setpoints, weather)
    radiation_flux_PAR_LampFlr = floor_PAR_absorbed_from_lamp(states, setpoints, weather)
    radiation_flux_NIR_LampFlr = floor_NIR_absorbed_from_lamp(states, setpoints, weather)
    radiation_flux_FIR_LampFlr = FIR_from_lamp_to_floor(states, setpoints, weather)

    radiation_flux_LampAir = lamp_radiation(states, setpoints, weather)
    sensible_heat_flux_LampCool = 0
    return (electrical_input_Lamp - radiation_flux_LampSky - radiation_flux_LampCov_in
                      - radiation_flux_LampThScr - radiation_flux_LampBlScr - sensible_heat_flux_LampAir
                      - radiation_flux_PAR_LampCanopy - radiation_flux_NIR_LampCanopy - radiation_flux_FIR_LampCanopy
                      - radiation_flux_LampPipe - radiation_flux_PAR_LampFlr - radiation_flux_NIR_LampFlr
                      - radiation_flux_FIR_LampFlr - radiation_flux_LampAir - sensible_heat_flux_LampCool)/cap_Lamp


def inter_lamps_temperature(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2 [2] [W m-2]
    cap_LampInt * intLamp_t = electrical_input_IntLampIn - sensible_heat_flux_IntLampAir - radiation_flux_PAR_IntLampCanopy
                            - radiation_flux_NIR_IntLampCanopy - radiation_flux_FIR_IntLampCanopy
    Returns:
    """
    cap_LampInt = 0
    electrical_input_IntLampIn = inter_lamp_electrical_input(states, setpoints)
    sensible_heat_flux_IntLampAir = sensible_heat_flux_between_inter_lamp_and_greenhouse_air(states)
    radiation_flux_PAR_IntLampCanopy = canopy_PAR_absorbed_from_inter_lamp(setpoints)
    radiation_flux_NIR_IntLampCanopy = canopy_NIR_absorbed_from_inter_lamp(setpoints)
    radiation_flux_FIR_IntLampCanopy = FIR_from_inter_lamp_to_canopy(states)
    return (electrical_input_IntLampIn - sensible_heat_flux_IntLampAir - radiation_flux_PAR_IntLampCanopy
            - radiation_flux_NIR_IntLampCanopy - radiation_flux_FIR_IntLampCanopy)/cap_LampInt


def greenhouse_air_vapor_pressure(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.10 / 8.10
    cap_vapor_Air * air_vapor_pressure = mass_vapor_flux_CanopyAir + mass_vapor_flux_FogAir
                                    + mass_vapor_flux_BlowAir − mass_vapor_flux_AirThScr − mass_vapor_flux_AirTop
                                    − mass_vapor_flux_AirOut − mass_vapor_flux_AirMech
    :return: The greenhouse air vapor pressure
    """
    cap_vapor_Air = air_compartment_water_vapor_capacity(states)
    mass_vapor_flux_CanopyAir = canopy_transpiration(states, setpoints, weather)
    mass_vapor_flux_FogAir = fogging_system_to_greenhouse_air_latent_vapor_flux(setpoints)
    mass_vapor_flux_BlowAir = heat_blower_to_greenhouse_air_vapor_flux(setpoints)
    mass_vapor_flux_AirThScr = greenhouse_air_to_thermal_screen_vapor_flux(setpoints, states)
    mass_vapor_flux_AirTop = greenhouse_air_to_above_thermal_screen_vapor_flux(states, setpoints, weather)
    mass_vapor_flux_AirOut = greenhouse_air_to_outdoor_vapor_flux(states, setpoints, weather)
    mass_vapor_flux_AirMech = greenhouse_air_to_mechanical_cooling_vapor_flux(states, setpoints)

    return (mass_vapor_flux_CanopyAir + mass_vapor_flux_FogAir
            + mass_vapor_flux_BlowAir - mass_vapor_flux_AirThScr - mass_vapor_flux_AirTop
            - mass_vapor_flux_AirOut - mass_vapor_flux_AirMech) / cap_vapor_Air


def top_compartment_vapor_pressure(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.11 / 8.11
    cap_vapor_Top * top_vapor_pressure = mass_vapor_flux_AirTop − mass_vapor_flux_TopCov_in − mass_vapor_flux_TopOut
    :return: The above thermal screen air vapor pressure
    """
    cap_vapor_Top = air_compartment_water_vapor_capacity(states)
    mass_vapor_flux_AirTop = greenhouse_air_to_above_thermal_screen_vapor_flux(states, setpoints, weather)
    mass_vapor_flux_TopCov_in = above_thermal_screen_to_internal_cover_vapor_flux(states)
    mass_vapor_flux_TopOut = above_thermal_screen_to_outdoor_vapor_flux(states, setpoints, weather)
    return (mass_vapor_flux_AirTop - mass_vapor_flux_TopCov_in - mass_vapor_flux_TopOut) / cap_vapor_Top


def greenhouse_air_co2(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.12 / 8.12
    cap_CO2_Air * air_CO2 = mass_CO2_flux_BlowAir + mass_CO2_flux_ExtAir
                          - mass_CO2_flux_AirCanopy - mass_CO2_flux_AirTop - mass_CO2_flux_AirOut
    :return: The greenhouse air CO2 concentration
    """
    cap_CO2_Air = coefs.Construction.air_height
    mass_CO2_flux_BlowAir = heat_blower_to_greenhouse_air_CO2_flux(setpoints)
    mass_CO2_flux_ExtAir = external_CO2_added(setpoints)
    mass_CO2_flux_AirCanopy = states.mass_CO2_flux_AirCanopy
    mass_CO2_flux_AirTop = greenhouse_air_and_above_thermal_screen_CO2_flux(states, setpoints, weather)
    mass_CO2_flux_AirOut = greenhouse_air_and_outdoor_CO2_flux(states, setpoints, weather)
    return (mass_CO2_flux_BlowAir + mass_CO2_flux_ExtAir
            - mass_CO2_flux_AirCanopy - mass_CO2_flux_AirTop - mass_CO2_flux_AirOut) / cap_CO2_Air


def top_compartment_air_CO2(setpoints: Setpoints, states: States, weather: Weather):
    """
    Equation 2.13 / 8.13
    cap_CO2_Top * top_CO2 = mass_CO2_flux_AirTop - mass_CO2_flux_TopOut
    :return: The above thermal screen air CO2 concentration
    """
    cap_CO2_Top = coefs.Construction.greenhouse_height - coefs.Construction.air_height
    mass_CO2_flux_AirTop = greenhouse_air_and_above_thermal_screen_CO2_flux(states, setpoints, weather)
    mass_CO2_flux_TopOut = above_thermal_screen_and_outdoor_CO2_flux(states, setpoints, weather)
    return (mass_CO2_flux_AirTop - mass_CO2_flux_TopOut) / cap_CO2_Top

