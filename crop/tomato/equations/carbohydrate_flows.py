from crop.tomato.equations.electron_transport import *
from crop.tomato.equations.inhibitions import *
from crop.tomato.equations.utils import *


def net_photosynthesis_rate():
    """
    Equation 9.10
    carbohydrate_mass_flow_AirBuf = M_CH2O * carbohydrates_saturation_photosynthesis_rate_inhibition
                                  * (gross_canopy_photosynthesis_rate - photosynthesis_process_photorespiration)
    Returns: The photosynthesis rate [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition()
    stomata_CO2_concentration = 0
    CO2_compensation_point = 0
    gross_canopy_photosynthesis_rate = canopy_level_photosynthesis_rate(stomata_CO2_concentration,
                                                                        CO2_compensation_point)
    photosynthesis_process_photorespiration_rate = photorespiration_rate(gross_canopy_photosynthesis_rate,
                                                                         stomata_CO2_concentration,
                                                                         CO2_compensation_point)
    return M_CH2O * carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * (gross_canopy_photosynthesis_rate - photosynthesis_process_photorespiration_rate)


def canopy_level_photosynthesis_rate(stomata_CO2_concentration, CO2_compensation_point):
    """
    Equations 9.12
    canopy_level_photosynthesis_rate = electron_transport_rate * (stomata_CO2_concentration - CO2_compensation_point) / (4 * (stomata_CO2_concentration + 2*CO2_compensation_point))
    Returns: gross canopy photosynthesis rate [µmol {CO2} m^-2 s^-1]
    """
    electron_transport_rate = electron_transport()
    return electron_transport_rate * (stomata_CO2_concentration - CO2_compensation_point) \
           / (4 * (stomata_CO2_concentration + 2 * CO2_compensation_point))


def photorespiration_rate(gross_canopy_photosynthesis_rate, stomata_CO2_concentration, CO2_compensation_point):
    """
    Equations 9.13
    photorespiration_rate = gross_canopy_photosynthesis_rate * CO2_compensation_point / stomata_CO2_concentration
    Args:
        gross_canopy_photosynthesis_rate:
        stomata_CO2_concentration:
        CO2_compensation_point:
    Returns: photorespiration rate [µmol {CO2} m^-2 s^-1]
    """
    return gross_canopy_photosynthesis_rate * CO2_compensation_point / stomata_CO2_concentration


def carbohydrate_flow_from_buffer_to_fruits():
    """
    Equations 9.24
    carbohydrate_flow_BufFruits = carbohydrates_saturation_photosynthesis_rate_inhibition_rate
                                * non_optimal_instantaneous_temperature_inhibition_rate
                                * non_optimal_24_hour_canopy_temperatures_inhibition_rate
                                * crop_development_stage_inhibition_rate
                                * growth_rate
                                * POTENTIAL_FRUIT_GROWTH_RATE_COEF
    Returns: carbohydrate flow from buffer to fruits [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition()
    non_optimal_instantaneous_temperature_inhibition_rate = non_optimal_instantaneous_temperature_inhibition()
    non_optimal_24_hour_canopy_temperatures_inhibition_rate = non_optimal_24_hour_canopy_temperatures_inhibition()
    crop_development_stage_inhibition_rate = crop_development_stage_inhibition()
    growth_rate = growth_rate_dependency_to_temperature()
    return carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * non_optimal_instantaneous_temperature_inhibition_rate \
           * non_optimal_24_hour_canopy_temperatures_inhibition_rate \
           * crop_development_stage_inhibition_rate \
           * growth_rate \
           * POTENTIAL_FRUIT_GROWTH_RATE_COEF


def carbohydrate_flow_from_buffer_to_leaves():
    """
    Equations 9.25
    carbohydrate_flow_BufLeaves = carbohydrates_saturation_photosynthesis_rate_inhibition_rate
                                * non_optimal_24_hour_canopy_temperatures_inhibition_rate
                                * growth_rate * POTENTIAL_LEAF_GROWTH_RATE_COEF
    Returns: carbohydrate flow from buffer to leaves [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition()
    non_optimal_24_hour_canopy_temperatures_inhibition_rate = non_optimal_24_hour_canopy_temperatures_inhibition()
    growth_rate = growth_rate_dependency_to_temperature()
    return carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * non_optimal_24_hour_canopy_temperatures_inhibition_rate \
           * growth_rate \
           * POTENTIAL_LEAF_GROWTH_RATE_COEF


def carbohydrate_flow_from_buffer_to_stem():
    """
    Equations 9.25
    carbohydrate_flow_BufStem = carbohydrates_saturation_photosynthesis_rate_inhibition_rate
                                * non_optimal_24_hour_canopy_temperatures_inhibition_rate
                                * growth_rate
                                * POTENTIAL_STEM_GROWTH_RATE_COEF
    Returns: carbohydrate flow from buffer to stem [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition()
    non_optimal_24_hour_canopy_temperatures_inhibition_rate = non_optimal_24_hour_canopy_temperatures_inhibition()
    growth_rate = growth_rate_dependency_to_temperature()
    return carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * non_optimal_24_hour_canopy_temperatures_inhibition_rate \
           * growth_rate \
           * POTENTIAL_STEM_GROWTH_RATE_COEF


def carbohydrates_flow_through_fruit_stages(jth: int):
    """
    Equations 9.34
    carbohydrates_flow_Fruit_j_Fruit_jplus = fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * carbohydrate_amount_Fruit_j
    Args:
        jth: stage number

    Returns: carbohydrates flow from fruit stage j to fruit stage j+1 [mg m^-2 s^-1]
    """
    carbohydrate_amount_Fruit_j = 0
    fruit_development_rate = fruit_development()
    return fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * carbohydrate_amount_Fruit_j


def carbohydrates_flow_from_buffer_to_first_fruit_stage():
    """
    Equations 9.35
    carbohydrates_flow_BuffFruit_1 = W_Pot_Fruit_1 * number_flow_BufFruit_1
    Returns: carbohydrates flow from buffer to fruit stage 1 [mg m^-2 s^-1]
    """
    pass


def carbohydrates_flow_from_buffer_to_fruit_stages(jth: int):
    """
    Equations 9.36
    carbohydrates_flow_BuffFruit_j = sum_carbohydrates_flow_BuffFruit_conversion_factor * number_fruit_j
                                    * fruit_growth_rate * (carbohydrates_flow_BuffFruit - carbohydrates_flow_BuffFruit_1)
    Returns: carbohydrates flow from buffer to fruit stage 1 [mg m^-2 s^-1]
    """
    sum_carbohydrates_flow_BufFruit_conversion_factor = sum_carbohydrates_flow_BufFruit_conversion()
    number_fruit_j = number_fruits[jth]
    fruit_growth_rate_j = fruit_growth(jth)
    carbohydrates_flow_BuffFruit = 0
    carbohydrates_flow_BuffFruit_1 = carbohydrates_flow_from_buffer_to_first_fruit_stage()
    return sum_carbohydrates_flow_BufFruit_conversion_factor * number_fruit_j \
           * fruit_growth_rate_j * (carbohydrates_flow_BuffFruit - carbohydrates_flow_BuffFruit_1)


def carbohydrates_flow_from_growth_respiration():
    """
    Equations 9.43
    carbohydrates_flow_BufAir =
    Returns: carbohydrates flow from growth respiration [mg m^-2 s^-1]
    """
    return FRUIT_GROWTH_RESPIRATION_COEF * carbohydrates_flow_BufFruit \
           + LEAF_GROWTH_RESPIRATION_COEF * carbohydrates_flow_BufLeaf \
           + STEM_GROWTH_RESPIRATION_COEF * carbohydrates_flow_BufStem


def carbohydrates_flow_from_fruit_maintenance_respiration():
    """
    Equations 9.45
    carbohydrates_flow_FruitAir = FRUIT_MAINTENANCE_RESPIRATION_COEF * carbohydrates_flow_BufFruit
    Returns: carbohydrates flow from fruit maintenance respiration [mg m^-2 s^-1]
    """
    return FRUIT_MAINTENANCE_RESPIRATION_COEF * carbohydrates_flow_BufFruit


def carbohydrates_flow_from_leaf_maintenance_respiration():
    """
    Equations 9.45
    carbohydrates_flow_LeafAir = LEAF_GROWTH_RESPIRATION_COEF * carbohydrates_flow_BufLeaf
    Returns: carbohydrates flow from leaf maintenance respiration [mg m^-2 s^-1]
    """
    return LEAF_MAINTENANCE_RESPIRATION_COEF * carbohydrates_flow_BufLeaf


def carbohydrates_flow_from_stem_maintenance_respiration():
    """
    Equations 9.45
    carbohydrates_flow_StemAir = STEM_GROWTH_RESPIRATION_COEF * carbohydrates_flow_BufStem
    Returns: carbohydrates flow from stem maintenance respiration [mg m^-2 s^-1]
    """
    return STEM_MAINTENANCE_RESPIRATION_COEF * carbohydrates_flow_BufStem


def leaf_harvest_rate():
    """
    Equations 9.47
    carbohydrates_flow_LeafHar =
    Returns: leaf harvest rate [mg m^-2 s^-1]
    """
    pass
