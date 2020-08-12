"""
These coefficients are imported from Table 8.1 in Vanthoor's thesis and from GreenLight code
NOTE: In the GreenLight model, there are no whitewash and shadow screen
"""


# TODO: import coefficients through a configuration file
class Coefficients(object):
    class Construction:
        ratio_GlobAir = 0.1  # The ratio of the global radiation which is absorbed by the greenhouse construction elements
        mean_greenhouse_cover_slope = 22  # Mean greenhouse cover slope
        cover_area = 9.0E4  # Surface of the cover including side-walls
        floor_area = 7.8E4  # Surface of the greenhouse floor
        c_HECin = 1.86  # Convective heat exchange parameter between cover and outdoor air that depends on the greenhouse shape
        c_HECout_1 = 2.8  # Convective heat exchange variables between cover and outdoor air which depend on the greenhouse shape
        c_HECout_2 = 1.2  # Convective heat exchange variables between cover and outdoor air which depend on the greenhouse shape
        c_HECout_3 = 1  # Convective heat exchange variables between cover and outdoor air which depend on the greenhouse shape
        air_height = 4.7  # Height of the greenhouse compartment below the thermal screen
        elevation_height = 1470  # The altitude of the greenhouse
        greenhouse_height = 5.1  # Mean height of the greenhouse
        C_Gh_d = 0.65  # Ventilation discharge coefficient depends on greenhouse shape
        C_Gh_w = 0.09  # Ventilation global wind pressure coefficient depends on greenhouse shape
        leakage_coef = 1E-4  # Greenhouse leakage_coef coefficient
        # side_wall_roof_vent_distance The vertical distance between mid-points of side wall and roof ventilation openings
        vent_vertical_dimension = 0.97  # The vertical dimension of a single ventilation opening

    class Ventilation:
        # eta_ShScrC_d Parameter that determines the effect of the movable shading screen on the discharge coefficient
        # eta_ShScrC_w Parameter that determines the effect of the movable shading screen on the global wind pressure coefficient
        A_Roof = 7.8E3  # 0.1*floor_area
        A_Side = 0
        porosity_InsScr = 1  # The porosity of the insect screens
        A_Roof_A_Flr = 0.18  # The specific roof ventilation area
        A_Side_A_Flr = 0  # The side ventilation area

    class Roof:
        roof_FIR_emission_coefficient = 0.85  # The FIR emission coefficient of the roof
        roof_density = 2.6E3  # Density of the roof layer
        roof_NIR_reflection_coefficient = 0.13  # The NIR reflection coefficient of the roof
        roof_PAR_reflection_coefficient = 0.13  # The PAR reflection coefficient of the roof
        roof_FIR_reflection_coefficient = 0.15  # The FIR reflection coefficient of the roof
        roof_NIR_transmission_coefficient = 0.85  # The NIR transmission coefficient of the roof
        roof_PAR_transmission_coefficient = 0.85  # The PAR transmission coefficient of the roof
        roof_FIR_transmission_coefficient = 0  # The FIR transmission coefficient of the roof
        roof_heat_conductivity = 1.05  # Thermal heat conductivity of the roof
        c_p_Rf = 0.84E3  # The specific heat capacity of the roof layer
        roof_thickness = 4E-3  # Thickness of the roof layer

    class Shadowscreen:
        # No shadowscreen
        shScr_NIR_reflection_coefficient = 0
        shScr_PAR_reflection_coefficient = 0
        shScr_FIR_reflection_coefficient = 0
        shScr_NIR_transmission_coefficient = 1
        shScr_PAR_transmission_coefficient = 1
        shScr_FIR_transmission_coefficient = 1

    class Thermalscreen:
        thScr_FIR_emission_coefficient = 0.44  # The FIR emission coefficient of the thermal screen
        thScr_density = 0.2E3  # Density of the thermal screen
        thScr_NIR_reflection_coefficient = 0.7  # The NIR reflection coefficient of the thermal screen
        thScr_PAR_reflection_coefficient = 0.7  # The PAR reflection coefficient of the thermal screen
        thScr_FIR_reflection_coefficient = 0.45  # The FIR reflection coefficient of the thermal screen
        thScr_NIR_transmission_coefficient = 0.25  # The NIR transmission coefficient of the thermal screen
        thScr_PAR_transmission_coefficient = 0.25  # The PAR transmission coefficient of the thermal screen
        thScr_FIR_transmission_coefficient = 0.11  # FIR transmission coefficient of the thermal screen
        c_pThScr = 1.8E3  # Specific heat capacity of the thermal screen
        thScr_thickness = 0.35E-3  # Thickness of the thermal screen
        thScr_flux_coefficient = 0.25E-3  # The thermal screen flux coefficient

    class Floor:
        floor_FIR_emission_coefficient = 1  # FIR emission coefficient of the floor
        floor_density = 2300  # Density of the floor
        floor_NIR_reflection_coefficient = 0.5  # NIR reflection coefficient of the floor
        floor_PAR_reflection_coefficient = 0.65  # PAR reflection coefficient of the floor
        floor_heat_conductivity = 1.7  # Thermal heat conductivity of the floor
        c_pFlr = 0.88E3  # Specific heat capacity of the floor
        floor_thickness = 0.02  # Thickness of the greenhouse floor

    class Soil:
        rho_c_p_So = 1.73E6  # The volumetric heat capacity of the soil
        soil_heat_conductivity = 0.85  # Thermal heat conductivity of the soil layers.
        soil_thicknesses = [0.04, 0.08, 0.16, 0.32, 0.64]  # The thickness of five soil layers

    class Heating:
        pipe_FIR_emission_coefficient = 0.88  # FIR emission coefficient of the heating pipes
        phi_external_pipe = 51E-3  # External diameter of the heating pipe
        phi_internal_pipe = 47E-3  # Internal diameter of the heating pipe
        pipe_length = 1.25  # Length of the heating pipes per square meter greenhouse

    class ActiveClimateControl:
        # cap_Fog Capacity of the fogging system
        # ventForced_air_flow Air flow capacity of the forced ventilation system
        cap_extCO2 = 4.3E5  # Capacity of the external CO2 source
        # perf_MechCool_coef Coefficient of performance of the mechanical cooling system
        # HEC_PasAir The convective heat exchange coefficient between the passive heat storage facility and the greenhouse air temperature
        # heat_cap_Blow Heat capacity of the heat blowers
        # heat_cap_Boil Thermal heat capacity of the boiler
        # heat_cap_Geo Heat capacity of the geothermal heat source
        # heat_cap_Ind Heat capacity of the industrial heat source
        # ele_cap_MechCool Electrical capacity of the mechanical cooling system

    class Lamp:
        # No lamps
        heat_capacity_lamp = 100  # Maximum heat capacity of lamps  [J K^-1 m^-2]
        electrical_capacity_lamp = 110  # electrical capacity of the lamps [W m-2]
        A_Lamp = 0.03  # surface area of the lamps per area of greenhouse floor [m m-2]
        lamp_NIR_reflection_coef = 0  # The NIR reflection coefficient of the vertical layer of the lamps
        lamp_PAR_reflection_coef = 0  # The PAR reflection coefficient of the vertical layer of the lamps
        lamp_FIR_reflection_coef = 0  # The FIR reflection coefficient of the vertical layer of the lamps
        lamp_NIR_transmission_coef = 0.97  # The NIR transmission coefficient of the vertical layer of the lamps
        lamp_PAR_transmission_coef = 0.97  # The PAR transmission coefficient of the vertical layer of the lamps
        lamp_FIR_transmission_coef = 0.97  # FIR transmission coefficient of the vertical layer of the lamps
        lamp_electrical_input_PAR_conversion = 0.36  # the conversion rate from electrical input to PAR output of the lamp
        lamp_electrical_input_NIR_conversion = 0.22  # the conversion rate from electrical input to NIR output of the lamp
        lamp_photons_per_joule = 5  # the amount of photons per joule within the PAR output of the lamps, which depends on the lamps' spectral output
        top_lamp_emission = 0.1  # the emissivity of the lamps towards the top
        bottom_lamp_emission = 0.9  # the emissivity of the lamps towards the bottom
        lamp_cool_energy = 0  # the amount of energy exported from the lamps by active cooling and removed from the greenhouse, expressed as a fraction of the electrical input [-]
        c_HEC_LampAir = 0.09  # the heat exchange coefficient between the lamps and the surrounding air [W K-1 m-2]

    class Interlight:
        # No lamps
        electrical_capacity_inter_lamp = 0  # electrical capacity of the lamps [W m-2]
        A_Inter_lamp = 0  # surface area of the lamps per area of greenhouse floor [m2 m-2]
        inter_lamp_electrical_input_PAR_conversion = 0  # fraction of lamp input converted to PAR [J(PAR) J-1 {electricity}]
        inter_lamp_electrical_input_NIR_conversion = 0  # fraction of lamp input converted to NIR [J(NIR) J-1 {electricity}]
        inter_lamp_photons_per_joule = 0  # the amount of photons per joule within the PAR output of the lamps, which depends on the lamps' spectral output [micro-mol{PAR} J-1 {PAR}]
        inter_lamp_emission = 0  # the emissivity of the lamps towards the top [-]
        heat_inter_lamp_capacity = 0  # heat capacity of lamps [J K^-1 m^-2]
        c_HEC_InterLampAir = 0  # the heat exchange coefficient between the lamps and the surrounding air [W K-1 m-2]

    class Blackoutscreen:
        blScr_PAR_transmission_coef = 0.01
        blScr_PAR_reflection_coef = 0.35
        blScr_FIR_emission_coef = 0.67  # FIR emissions coefficient of the blackout screen

    class GrowPipe:
        phi_external_pipe = 0.035  # External diameter of the grow pipes [m]
        phi_internal_pipe = 0.0338  # Internal diameter of the grow pipes [m]
        pipe_length = 1.655  # Length of the grow pipes per square meter greenhouse floor [m m-2]
        groPipe_FIR_emission_coef = 0  # FIR emissions coefficient of the blackout screen
