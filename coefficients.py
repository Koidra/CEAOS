from math import inf

"""
NOTE: In the GreenLight model, there are no whitewash and shadow screen
"""


class Coefficients(object):
    def __init__(self):
        pass

    class Outside:
        canopy_air_convective_heat_exchange_coefficient = 5  # Convective heat exchange coefficient from the canopy leaf to the greenhouse air
        evaporation_latent_heat = 2.45E6  # Latent heat of evaporation
        sigma = 5.670E-8  # Stefan Boltzmann constant
        can_FIR_emission_coefficient = 1  # FIR emission coefficient of the canopy
        sky_FIR_emission_coefficient = 1  # FIR emission coefficient of the sky
        ratio_GlobNIR = 0.5  # Ratio between NIR and the outside global radiation
        ratio_GlobPAR = 0.5  # Ratio between PAR and the outside global radiation
        eta_HeatCO2 = 0.057  # Amount of CO2 which is released when 1 Joule of sensible energy is produced by the heat blower
        eta_HeatVap = 4.43E-8  # Amount of vapor which is released when 1 Joule of sensible energy is produced by the heat blower
        eta_mg_ppm = 0.554  # CO2 conversion factor from mg m -3 to ppm.
        eta_Roof_Thr = 0.9  # The ratio between the roof vent area and total ventilation area above no chimney effect was assumed
        density_Air0 = 1.20  # Density of the air at sea level
        canopy_PAR_reflection_coefficient = 0.07  # The PAR reflection coefficient
        canopy_NIR_reflection_coefficient = 0.35  # The NIR reflection coefficient of the top of the canopy
        steel_density = 7850  # Density of steel
        water_density = 1E3  # Density of water
        gamma = 65.8  # Psychrometric constant
        omega = 1.99E-7  # The yearly frequency to calculate the soil temperature
        cap_Leaf = 1.2E3  # Heat capacity of a square meter canopy leaves
        c_evap1 = 4.30  # Coefficient of the stomatal resistance model to account for radiation effect
        c_evap2 = 0.54  # Coefficient of the stomatal resistance model to account for radiation effect
        c_day_evap3 = 6.1E-7  # Coefficient of the stomatal resistance model to account CO2 effect
        c_night_evap3 = 1.1E-11  # Coefficient of the stomatal resistance model to account CO2 effect
        c_day_evap4 = 4.3E-6  # Coefficient of the stomatal resistance model to account for vapor pressure difference
        c_night_evap4 = 5.2E-6  # Coefficient of the stomatal resistance model to account for vapor pressure difference
        c_pAir = 1E3  # Specific heat capacity of the air
        c_pSteel = 0.64E3  # Specific heat capacity of steel
        c_pWater = 4.18E3  # Specific heat capacity of water
        g = 9.81  # The acceleration of gravity
        canopy_PAR_extinction_coefficient = 0.7  # PAR extinction coefficient of the canopy
        floor_PAR_extinction_coefficient = 0.7  # PAR extinction coefficient of the canopy when PAR is reflected from the floor
        canopy_NIR_extinction_coefficient = 0.27  # Extinction coefficient of the canopy for NIR
        canopy_FIR_extinction_coefficient = 0.94  # Extinction coefficient of the canopy for FIR
        M_Air = 28.96  # Molar mass of air
        M_Water = 18  # Molar mass of water
        R = 8.314E3  # Molar gas constant
        R_Can_SP = 5  # The radiation value above the canopy when the night becomes day and vice versa
        r_b = 275  # Boundary layer resistance of the canopy for vapor transport
        r_s_min = 82.0  # The minimum canopy resistance for transpiration
        s_r_s = -1  # The slope of the differentiable switch for the stomatical resistance model
        s_MV12 = -0.1  # The slope of the differentiable switch function for vapor pressure differences

    class Greenhouse:
        class Construction:
            eta_GlobAir = 0.1  # The ratio of the global radiation which is absorbed by the greenhouse construction elements
            mean_greenhouse_cover_slope = 22  # Mean greenhouse cover slope
            cover_surface = 9.0E4  # Surface of the cover including side-walls
            floor_surface = 7.8E4  # Surface of the greenhouse floor
            c_HECin = 1.86  # Convective heat exchange parameter between cover and outdoor air that depends on the greenhouse shape
            c_HECout_1 = 2.8  # Convective heat exchange variables between cover and outdoor air which depend on the greenhouse shape
            c_HECout_2 = 1.2  # Convective heat exchange variables between cover and outdoor air which depend on the greenhouse shape
            c_HECout_3 = 1  # Convective heat exchange variables between cover and outdoor air which depend on the greenhouse shape
            air_height = 4.7  # Height of the greenhouse compartment below the thermal screen
            elevation_height = 1470  # The altitude of the greenhouse
            greenhouse_height = 5.1  # Mean height of the greenhouse
            C_Gh_d = 0.65  # Ventilation discharge coefficient depends on greenhouse shape
            C_Gh_w = 0.09  # Ventilation global wind pressure coefficient depends on greenhouse shape
            leakage = 1E-4  # Greenhouse leakage coefficient
            # h_SideRoof The vertical distance between mid-points of side wall and roof ventilation openings
            h_Vent = 0.97  # The vertical dimension of a single ventilation opening

        class Ventilation:
            # eta_ShScrC_d Parameter that determines the effect of the movable shading screen on the discharge coefficient
            # eta_ShScrC_w Parameter that determines the effect of the movable shading screen on the global wind pressure coefficient
            A_Roof = 7.8E3  # 0.1*floor_surface
            A_Side = 0
            sigma_InsScr = 1  # The porosity of the insect screens
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

        # TODO: What is the differences between Shadowscreen and Blackoutscreen
        class Blackoutscreen:
            blScr_PAR_transmission_coefficient = 0.01
            blScr_PAR_reflection_coefficient = 0.35

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
            # phi_Fog Capacity of the fogging system
            # phi_VentForced Air flow capacity of the forced ventilation system
            extCO2_capacity = 4.3E5  # Capacity of the external CO2 source
            # COP_MechCool Coefficient of performance of the mechanical cooling system
            # HEC_PasAir The convective heat exchange coefficient between the passive heat storage facility and the greenhouse air temperature
            # P_Blow Heat capacity of the heat blowers
            # P_Boil Thermal heat capacity of the boiler
            # P_Geo Heat capacity of the geothermal heat source
            # P_Ind Heat capacity of the industrial heat source
            # P_MechCool Electrical capacity of the mechanical cooling system

        class Lamp:
            # No lamps
            thetaLampMax = 0  # Maximum intensity of lamps
            eta_LampPAR = 0  # fraction of lamp input converted to PAR
            eta_LampNIR = 0  # fraction of lamp input converted to NIR

        class Interlight:
            # No lamps
            thetaIntLampMax = 0  # Maximum intensity of lamps
            eta_IntLampPAR = 0  # fraction of lamp input converted to PAR
            eta_IntLampNIR = 0  # fraction of lamp input converted to NIR
