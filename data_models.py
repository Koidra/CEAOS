from typing import NamedTuple


class Setpoints(NamedTuple):
    U_Blow: float  # Heat blower control
    U_Boil: float  # Boiler control
    U_MechCool: float  # Control of the mechanical cooling
    U_Fog: float  # Control of fogging system
    U_Roof: float  # Control of the roof ventilators
    U_Side: float  # Control of the side ventilators
    U_VentForced: float  # Control of the forced ventilation
    U_ExtCO2: float  # Control of the CO2-input from an external source
    U_ShScr: float  # Control of the external shading screen
    U_ThScr: float  # Control of the thermal screen
    U_Ind: float  # Control of the heat input from industry
    U_Geo: float  # Control of the heat input from geothermal source
    U_Lamp: float  # the switching of the top-lights
    U_IntLamp: float  # the switching of the inter-lights
    U_BoilGro: float  # the valve opening between the boiler and the grow pipes
    U_BlScr: float  # the opening of the blackout screen


class States(NamedTuple):
    pipe_t: float  # pipe temperature
    canopy_t: float  # canopy temperature
    air_t: float  # greenhouse air temperature
    internal_cov_t: float  # internal cover temperature
    external_cov_t: float  # external cover temperature
    thermal_screen_t: float  # thermal screen temperature
    above_thermal_screen_t: float  # above thermal screen temperature a.k.a top compartment temperature
    floor_t: float  # floor temperature
    soil_j_t: [float, float, float, float, float]  # soil layer temperatures
    blScr_t: float
    groPipe_t: float
    lamp_t: float
    intLamp_t: float
    air_CO2: float  # CO2 in greenhouse air
    above_thermal_screen_CO2: float  # CO2 in top compartment air
    air_vapor_pressure: float
    top_vapor_pressure: float
    # Recheck these vars
    leaf_area_index: float
    mechcool_t: float  # Mechanical cooling system temperature
    mass_CO2_flux_AirCanopy: float  # CO2 flux from greenhouse air to canopy


class Weather(NamedTuple):
    outdoor_global_rad: float  # outdoor global radiation
    outdoor_t: float  # outdoor temperature
    sky_t: float  # sky temperature
    soil_out_t: float  # deep out soil temperature
    outdoor_CO2: float  # outdoor CO2
    outdoor_vp: float  # outdoor vapor pressure
    v_Wind: float  # wind velocity
