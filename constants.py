"""
These constants are imported from Table 8.1 and 9.1 in Vanthoor's thesis
"""

# Convective heat exchange coefficient from the canopy leaf to the greenhouse air.
# Unit: W m^-2 K^-1.
# Ref: De Zwart (1996)
CANOPY_AIR_CONVECTIVE_HEAT_EXCHANGE_COEF = 5

# Latent heat of evaporation.
# Unit: J kg^-1 {water}
EVAPORATION_LATENT_HEAT = 2.45E6

# Stefan Boltzmann constant.
# Unit: W m^-2 K^-4
BOLTZMANN = 5.670E-8

# FIR emission coefficient of the canopy.
# Unit: [-].
# Ref: Stanghellini (1987): leaves are considered as black bodies
CANOPY_FIR_EMISSION_COEF = 1

# FIR emission coefficient of the sky.
# Unit: [-]. Ref: by definition
SKY_FIR_EMISSION_COEF = 1

# Ratio between NIR and the outside global radiation.
# Unit: [-].
# Ref: Monteith (1973): UV is attributed to NIR.
RATIO_GLOBALNIR = 0.5

# Ratio between PAR and the outside global radiation.
# Unit: [-].
# Ref: Monteith (1973)
RATIO_GLOBALPAR = 0.5

# Amount of CO2 which is released when 1 Joule of sensible energy is produced by the heat blower.
# Unit: mg {CO2} J^-1
ETA_HEATCO2 = 0.057

# Amount of vapor which is released when 1 Joule of sensible energy is produced by the heat blower.
# Unit: kg {vapour} J-1
ETA_HEATVAP = 4.43E-8

# CO2 conversion factor from mg m -3 to ppm.
# Unit: ppm mg^-1 m^3
ETA_MG_PPM = 0.554

# The ratio between the roof vent area and total ventilation area above no chimney effect was assumed.
# Unit: [-].
# Ref: Assumed by Vanthoor
ETA_ROOF_THR = 0.9

# Density of the air at sea level.
# Unit: kg m^-3
DENSITY_AIR0 = 1.20

# The PAR reflection coefficient of the top of the canopy.
# Unit: [-].
# Ref: Marcelis et al. (1998)
CANOPY_PAR_REFLECTION_COEF = 0.07

# The NIR reflection coefficient of the top of the canopy.
# Unit: [-].
# Ref: De Zwart (1996): based on absorption of diffuse NIR
CANOPY_NIR_REFLECTION_COEF = 0.35

# Density of steel.
# Unit: kg m^-3
STEEL_DENSITY = 7850

# Density of water.
# Unit: kg m^-3
WATER_DENSITY = 1E3

# Psychrometric constant.
# Unit: Pa K^-1
GAMMA = 65.8

# The yearly frequency to calculate the soil temperature.
# Unit: s^-1
OMEGA = 1.99E-7

# Heat capacity of a square meter canopy leaves.
# Unit: J K^-1 m^-2 {leaf}.
# Ref: Stanghellini (1987)
CAP_LEAF = 1.2E3

# Coefficient of the stomatal resistance model to account for radiation effect.
# Unit: W m^-2.
# Ref: Stanghellini (1987)
C_EVAP1 = 4.30

# Coefficient of the stomatal resistance model to account for radiation effect.
# Unit: W m^-2.
# Ref: Stanghellini (1987)
C_EVAP2 = 0.54

# Coefficient of the stomatal resistance model to account CO2 effect.
# Unit: ppm^-2.
# Ref: Stanghellini (1987)
C_DAY_EVAP3 = 6.1E-7

# Coefficient of the stomatal resistance model to account CO2 effect.
# Unit: ppm^-2.
# Ref: Stanghellini (1987)
C_NIGHT_EVAP3 = 1.1E-11

# Coefficient of the stomatal resistance model to account for vapor pressure difference.
# Unit: Pa^-2.
# Ref: Stanghellini (1987)
C_DAY_EVAP4 = 4.3E-6

# Coefficient of the stomatal resistance model to account for vapor pressure difference.
# Unit: Pa^-2.
# Ref: Stanghellini (1987)
C_NIGHT_EVAP4 = 5.2E-6

# Specific heat capacity of the air.
# Unit: J K^-1 kg^-1
C_PAIR = 1E3

# Specific heat capacity of steel.
# Unit: J K^-1 kg^-1
C_PSTEEL = 0.64E3

# Specific heat capacity of water.
# Unit: J K^-1 kg^-1
C_PWATER = 4.18E3

# The acceleration of gravity. Unit: m s^-2
GRAVITY = 9.81

# PAR extinction coefficient of the canopy.
# Unit: [-].
# Ref: Marcelis et al. (1998)
CANOPY_PAR_EXTINCTION_COEF = 0.7

# PAR extinction coefficient of the canopy when PAR is reflected from the floor.
# Unit: [-].
# Ref: Assumed for white mulching by Vanthoor
FLOOR_PAR_EXTINCTION_COEF = 0.7

# Extinction coefficient of the canopy for NIR.
# Unit: [-].
# Ref: Based on absorption of diffuse NIR of De Zwart (1996)
CANOPY_NIR_EXTINCTION_COEF = 0.27

# Extinction coefficient of the canopy for FIR.
# Unit: [-].
# Ref: De Zwart (1996)
CANOPY_FIR_EXTINCTION_COEF = 0.94

# Molar mass of air.
# Unit: kg kmol^-1
M_AIR = 28.96

# Molar mass of water.
# Unit: kg kmol^-1
M_WATER = 18

# Molar gas constant.
# Unit: J kmol^-1 K^-1
M_GAS = 8.314E3

# The radiation value above the canopy when the night becomes day and vice versa.
# Unit: W m^-2
RAD_CANOPY_SETPOINT = 5

# Boundary layer resistance of the canopy for vapor transport.
# Unit: s m^-1.
# Ref: Mean value based on Stanghellini (Stanghellini, 1987)
BOUNDARY_LAYER_RESISTANCE = 275

# The minimum canopy resistance for transpiration.
# Unit: s m^-1.
# Ref: Stanghellini (1987)
MIN_CANOPY_TRANSPIRATION_RESISTANCE = 82.0

# The slope of the differentiable switch for the stomatical resistance model.
# Unit: m W^-2
S_R_S = -1

# The slope of the differentiable switch function for vapor pressure differences.
# Unit: Pa^-1
S_MV12 = -0.1

# Conversion factor of sun's PAR from J to
# Unit: umol{photons} J^{-1}
# Ref: GreenLight
PAR_J_TO_UMOL_SUN_RATE = 4.6