from typing import NamedTuple


class Setpoints(NamedTuple):
    U_Blow: float  # Heat blower control
    U_Boil: float  # Boiler control
    U_MechCool: float  # Control of the mechanical cooling
    U_Fog: float  # Control of fogging system
    U_Roof: float  # Control of the roof ventilators
    U_Side: float  # Control of the side ventilators
    U_VentForced: float  # Control of the forced ventilation
    U_Extco2: float  # Control of the CO2-input from an external source
    U_ShScr: float  # Control of the external shading screen
    U_ThScr: float  # Control of the thermal screen
    U_Ind: float  # Control of the heat input from industry
    U_Geo: float  # Control of the heat input from geothermal source
    U_Lamp: float  # the switching of the top-lights
    U_IntLamp: float  # the switching of the inter-lights
    U_BoilGro: float  # the valve opening between the boiler and the grow pipes
    U_BlScr: float  # the opening of the blackout screen


class ClimateStates(NamedTuple):
    t_Pipe: float  # pipe temperature
    t_Canopy: float  # canopy temperature
    t_Air: float  # greenhouse air temperature
    t_Cov_internal: float  # internal cover temperature
    t_Cov_external: float  # external cover temperature
    t_ThScr: float  # thermal screen temperature
    t_AboveThScr: float  # above thermal screen temperature a.k.a top compartment temperature
    t_Floor: float  # floor temperature
    t_Soil: [float, float, float, float, float]  # soil layer temperatures
    t_BlScr: float  # blackout screen temperatures
    t_GrowPipe: float  # grow pipe temperatures
    t_Lamp: float
    t_IntLamp: float
    co2_Air: float  # CO2 in greenhouse air
    co2_AboveThScr: float  # CO2 in top compartment air
    vapor_pressure_Air: float
    vapor_pressure_AboveThScr: float
    # Recheck these vars
    leaf_area_index: float
    t_MechCool: float  # Mechanical cooling system temperature
    mass_co2_flux_AirCanopy: float  # CO2 flux from greenhouse air to canopy
    PAR_Canopy: float


class Weather(NamedTuple):
    outdoor_global_rad: float  # outdoor global radiation
    t_Outdoor: float  # outdoor temperature
    t_Sky: float  # sky temperature
    t_Soil_Out: float  # deep out soil temperature
    co2_outdoor: float  # outdoor CO2
    vapor_pressure_outdoor: float  # outdoor vapor pressure
    v_Wind: float  # wind velocity
