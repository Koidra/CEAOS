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

# The conversion factor from photons to electrons including an efficiency term.
# Unit: µmol {e-} µmol^-1 {photons}.
# Ref: Farquhar et al. (1980)
PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR = 0.385

# Degree of curvature of the electron transport rate.
# Unit: [-].
# Ref: Farquhar et al. (1980)
ELECTRON_TRANSPORT_RATE_CURVATURE = 0.7

# Conversion factor from carbohydrate to dry matter.
# Unit: mg {DM} mg^-1 {CH2O}.
# Ref: No lignification assumed
CARBOHYDRATE_TO_DRY_MATTER_CONVERSION = 1

# Conversion factor from the greenhouse air CO2-concentration to the CO2-concentration in the stomata.
# Unit: µmol {CO2} mol^-1 {air}.
# Ref: Evans and Farquhar (1991)
GREENHOUSE_AIR_TO_STOMATA_CO2_CONCENTRATION = 0.67

# Conversion factor from global radiation to PAR.
# Unit: µmol {photons} J^-1.
# Ref:  De Zwart (1996)
GLOBAL_RADIATION_TO_PAR_CONVERSION = 2.3

# The time constant to calculate the 24 hour mean temperature.
# Unit: s.
# Ref: See section 9.7.1
DAY_MEAN_TEMP_TIME_CONSTANT = 86400

# Light transmission of the greenhouse cover.
# Unit: [-].
# Ref: Measured for Dutch growers
COVER_LIGHT_TRANSMISSION = 0.78

# The effect of canopy temperature on the CO2 compensation point.
# Unit: µmol {CO2} mol^-1 {air} K^-1
# Ref: Farquhar (1988)
CANOPY_TEMPERATURE_EFFECT = 1.7

# Maximum fruit set regression coefficient 1.
# Unit: fruits plant^-1 s^-1
# Ref: De Koning (1994)
MAX_FRUIT_SET_REGRESSION_COEF_1 = -1.71e-7

# Maximum fruit set regression coefficient 2.
# Unit: fruits plant^-1 s^-1 °C^-1
# Ref: De Koning (1994)
MAX_FRUIT_SET_REGRESSION_COEF_2 = 7.31e-7

# Fruit development rate coefficient 1.
# Unit: s^-1
# Ref: De Koning (1994)
FRUIT_DEVELOPMENT_RATE_COEF_1 = -7.64e-9

# Fruit development rate coefficient 2.
# Unit: s^-1 °C^-1
# Ref: De Koning (1994)
FRUIT_DEVELOPMENT_RATE_COEF_2 = 1.16e-8

# Fruit maintenance respiration coefficient.
# Unit: mg {CH2O} mg^-1 {CH2O} s^-1
# Ref: Heuvelink (1996)
FRUIT_MAINTENANCE_RESPIRATION_COEF = 1.16e-7

# Fruit growth respiration coefficient.
# Unit: [-]
# Ref: Heuvelink (1996)
FRUIT_GROWTH_RESPIRATION_COEF = 0.27

# Leaf maintenance respiration coefficient.
# Unit: mg {CH2O} mg^-1 {CH2O} s^-1
# Ref: Heuvelink (1996)
LEAF_MAINTENANCE_RESPIRATION_COEF = 3.47e-7

# Leaf growth respiration coefficient.
# Unit: [-]
# Ref: Heuvelink (1996)
LEAF_GROWTH_RESPIRATION_COEF = 0.28

# Stem maintenance respiration coefficient.
# Unit: mg {CH2O} mg^-1 {CH2O} s^-1
# Ref: Heuvelink (1996)
STEM_MAINTENANCE_RESPIRATION_COEF = 1.47e-7

# Stem growth respiration coefficient.
# Unit: [-]
# Ref: Heuvelink (1996)
STEM_GROWTH_RESPIRATION_COEF = 0.30

# Regression coefficient in maintenance respiration function.
# Unit: s
# Ref: Heuvelink (1996)
MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF = 2.85e6

# Maximum buffer capacity.
# Unit: mg {CH2O} m^-2
# Ref: Assumed by Vanthoor
MAX_BUFFER_CAP = 20e3

# Minimum amount of carbohydrates in the buffer.
# Unit: mg {CH2O} m^-2
# Ref: Assumed by Vanthoor
MIN_CARBOHYDRATES_IN_BUFFER = 1e3

# Activation energy for JPOT calculation.
# Unit: J mol^-1
# Ref: Farquhar et al. (1980)
ACTIVATION_ENERGY_JPOT = 37e3

# Potential fruit dry weight at harvest.
# Unit: mg {CH2O} fruit^-1
# Ref: De Koning (1994)
POTENTIAL_FRUIT_DRY_WEIGHT = 1e4

# Deactivation energy for JPOT calculation.
# Unit: J mol^-1
# Ref: Farquhar et al. (1980)
DEACTIVATION_ENERGY_JPOT = 22e4

# Maximal rate of electron transport at 25°C for the leaf.
# Unit: µmol {e-} m^-2 {leaf} s^-1
# Ref: Farquhar et al. (1980)
MAX_LEAF_ELECTRON_TRANSPORT_RATE = 210

# The gain of the process to calculate the 24 hour mean temperature.
# Unit: [-]
# Ref: See section 9.7.1
PROCESS_GAIN = 1

# Molar mass of CH2O.
# Unit: mg {CH2O} µmol-1 {CH2O}
M_CH2O = 30e-3

# Plant density in the greenhouse.
# Unit: plants m^-2
# Ref: Measured for Dutch growers
PLANT_DENSITY = 2.5

# Number of fruit development stages.
# Unit: [-]
# Ref: Assumption
FRUIT_DEVELOPMENT_STAGES_NUM = 50

# Carbohydrate flow from buffer to the fruits above which fruit set is maximal.
# Unit: mg {CH2O} m^-2 s^-1
# Ref: Assumption
BUFFER_TO_FRUITS_CARBOHYDRATE = 0.1

# Potential fruit growth rate coefficient at 20°C.
# Unit: mg {CH2O} m^-2 s^-1
# Ref: See section 9.7.3
POTENTIAL_FRUIT_GROWTH_RATE_COEF = 0.328

# Potential leaf growth rate coefficient at 20°C.
# Unit: mg {CH2O} m^-2 s^-1
# Ref: See section 9.7.3
POTENTIAL_LEAF_GROWTH_RATE_COEF = 0.095

# Potential stem growth rate coefficient.
# Unit: mg {CH2O} m^-2 s^-1
# Ref: See section 9.7.3
POTENTIAL_STEM_GROWTH_RATE_COEF = 0.074

# Entropy term for JPOT calculation.
# Unit: J mol^-1 K^-1
# Ref: Farquhar et al. (1980)
ENTROPY_TERM_JPOT = 710

# Specific leaf area index.
# Unit: m^2 {leaf} mg^-1 {CH2O}
# Ref: Assumed to be constant, Heuvelink (1996)
SPECIFIC_LEAF_AREA_INDEX = 2.66e-5

# Reference temperature for JPOT calculation.
# Unit: K
REFERENCE_TEMPERATURE_JPOT = 298.15

# Temperature sum when fruit growth rate is at full potential.
# Unit: °C
# Ref: Based upon Eq. (9.32)
SUM_END_T = 1035

# Base temperature for 24 hour mean crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
BASE_24_T = 12

# First optimal temperature for 24 hour mean crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
OPT1_24_T = 18

# Second optimal temperature for 24 hour mean crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
OPT2_24_T = 22

# Maximum temperature for 24 hour mean crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
MAX_24_T = 27

# Base temperature for instantaneous crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
BASE_INSTANTANEOUS_T = 6

# First optimal temperature for instantaneous crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
OPT1_INSTANTANEOUS_T = 14

# Second optimal temperature for instantaneous crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
OPT2_INSTANTANEOUS_T = 28

# Maximum temperature for instantaneous crop growth inhibition.
# Unit: °C
# Ref: See section 9.7.3
MAX_INSTANTANEOUS_T = 40

# Q10 value of temperature effect on maintenance respiration.
# Unit: [-]
# Ref: Heuvelink (1996)
Q10_M = 2


# Leaves are pruned back to this LAI
# Unit: [m^2 m^-2]
# Ref: GreenLight
LAI_Max = 3