from .CO2_concentration import co2_concentration_inside_stomata, co2_compensation
from .electron_transport import electron_transport
from .fruit_flow import fruit_set_of_first_development_stage
from .inhibitions import *
from .utils import *


def net_photosynthesis_rate(carbohydrate_amount_Buf, carbohydrate_amount_Leaf, air_co2, canopy_t, PAR_Canopy):
    """
    Equation 9.10
    carbohydrate_mass_flow_AirBuf = M_CH2O * carbohydrates_saturation_photosynthesis_rate_inhibition
                                  * (gross_canopy_photosynthesis_rate - photosynthesis_process_photorespiration)
    Returns: The photosynthesis rate [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition(
        carbohydrate_amount_Buf)
    stomata_co2_concentration = co2_concentration_inside_stomata(air_co2)
    co2_compensation_point = co2_compensation(carbohydrate_amount_Leaf, canopy_t)
    gross_canopy_photosynthesis_rate = canopy_level_photosynthesis_rate(carbohydrate_amount_Leaf, canopy_t,
                                                                        stomata_co2_concentration,
                                                                        co2_compensation_point, PAR_Canopy)
    return M_CH2O * carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * (gross_canopy_photosynthesis_rate
              - photorespiration_rate(gross_canopy_photosynthesis_rate,
                                      stomata_co2_concentration,
                                      co2_compensation_point))


def canopy_level_photosynthesis_rate(carbohydrate_amount_Leaf, canopy_t, stomata_co2_concentration, co2_compensation_point, PAR_Canopy):
    """
    Equations 9.12
    canopy_level_photosynthesis_rate = electron_transport_rate * (stomata_CO2_concentration - CO2_compensation_point) / (4 * (stomata_CO2_concentration + 2*CO2_compensation_point))
    Returns: gross canopy photosynthesis rate [µmol {CO2} m^-2 s^-1]
    """
    electron_transport_rate = electron_transport(carbohydrate_amount_Leaf, canopy_t, PAR_Canopy)
    return electron_transport_rate * (stomata_co2_concentration - co2_compensation_point) \
           / (4 * (stomata_co2_concentration + 2 * co2_compensation_point))


def photorespiration_rate(gross_canopy_photosynthesis_rate, stomata_co2_concentration, co2_compensation_point):
    """
    Equations 9.13
    photorespiration_rate = gross_canopy_photosynthesis_rate * CO2_compensation_point / stomata_CO2_concentration
    Args:
        gross_canopy_photosynthesis_rate:
        stomata_co2_concentration:
        co2_compensation_point:
    Returns: photorespiration rate [µmol {CO2} m^-2 s^-1]
    """
    return gross_canopy_photosynthesis_rate * co2_compensation_point / stomata_co2_concentration


def carbohydrate_flow_from_buffer_to_fruits(carbohydrate_amount_Buf, canopy_t, sum_canopy_t, last_24_canopy_t):
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
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition(carbohydrate_amount_Buf)
    non_optimal_instantaneous_temperature_inhibition_rate = non_optimal_instantaneous_temperature_inhibition(canopy_t)
    non_optimal_24_hour_canopy_temperatures_inhibition_rate = non_optimal_24_hour_canopy_temperatures_inhibition(last_24_canopy_t)
    crop_development_stage_inhibition_rate = crop_development_stage_inhibition(sum_canopy_t)
    growth_rate = growth_rate_dependency_to_temperature(last_24_canopy_t)
    return carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * non_optimal_instantaneous_temperature_inhibition_rate \
           * non_optimal_24_hour_canopy_temperatures_inhibition_rate \
           * crop_development_stage_inhibition_rate \
           * growth_rate \
           * POTENTIAL_FRUIT_GROWTH_RATE_COEF


def carbohydrate_flow_from_buffer_to_leaves(carbohydrate_amount_Buf, last_24_canopy_t):
    """
    Equations 9.25
    carbohydrate_flow_BufLeaf = carbohydrates_saturation_photosynthesis_rate_inhibition_rate
                                * non_optimal_24_hour_canopy_temperatures_inhibition_rate
                                * growth_rate * POTENTIAL_LEAF_GROWTH_RATE_COEF
    Returns: carbohydrate flow from buffer to leaves [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition(carbohydrate_amount_Buf)
    non_optimal_24_hour_canopy_temperatures_inhibition_rate = non_optimal_24_hour_canopy_temperatures_inhibition(last_24_canopy_t)
    growth_rate = growth_rate_dependency_to_temperature(last_24_canopy_t)
    return carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * non_optimal_24_hour_canopy_temperatures_inhibition_rate \
           * growth_rate \
           * POTENTIAL_LEAF_GROWTH_RATE_COEF


def carbohydrate_flow_from_buffer_to_stem(carbohydrate_amount_Buf, last_24_canopy_t):
    """
    Equations 9.25
    carbohydrate_flow_BufStem = carbohydrates_saturation_photosynthesis_rate_inhibition_rate
                                * non_optimal_24_hour_canopy_temperatures_inhibition_rate
                                * growth_rate
                                * POTENTIAL_STEM_GROWTH_RATE_COEF
    Returns: carbohydrate flow from buffer to stem [mg m^-2 s^-1]
    """
    carbohydrates_saturation_photosynthesis_rate_inhibition_rate = carbohydrates_saturation_photosynthesis_rate_inhibition(carbohydrate_amount_Buf)
    non_optimal_24_hour_canopy_temperatures_inhibition_rate = non_optimal_24_hour_canopy_temperatures_inhibition(last_24_canopy_t)
    growth_rate = growth_rate_dependency_to_temperature(last_24_canopy_t)
    return carbohydrates_saturation_photosynthesis_rate_inhibition_rate \
           * non_optimal_24_hour_canopy_temperatures_inhibition_rate \
           * growth_rate \
           * POTENTIAL_STEM_GROWTH_RATE_COEF


def carbohydrate_flow_through_fruit_stages(jth: int, carbohydrate_amount_Fruits, last_24_canopy_t):
    """
    Equations 9.34
    carbohydrate_flow_Fruit_j_Fruit_jplus = fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * carbohydrate_amount_Fruit_j
    Args:
        jth: stage number

    Returns: carbohydrates flow from fruit stage j to fruit stage j+1 [mg m^-2 s^-1]
    """
    carbohydrate_amount_Fruit_j = carbohydrate_amount_Fruits[jth - 1]
    fruit_development_rate = fruit_development(last_24_canopy_t)
    return fruit_development_rate * FRUIT_DEVELOPMENT_STAGES_NUM * carbohydrate_amount_Fruit_j


def carbohydrate_flow_from_buffer_to_first_fruit_stage(carbohydrate_amount_Buf, canopy_t, sum_canopy_t, last_24_canopy_t):
    """
    Equations 9.35
    carbohydrate_flow_BufFruit_1 = W_Pot_Fruit_1 * number_flow_BufFruit_1
    Returns: carbohydrates flow from buffer to fruit stage 1 [mg m^-2 s^-1]
    """
    carbohydrate_flow_BufFruits = carbohydrate_flow_from_buffer_to_fruits(carbohydrate_amount_Buf, canopy_t,
                                                                          sum_canopy_t, last_24_canopy_t)
    number_flow_BufFruit_1 = fruit_set_of_first_development_stage(last_24_canopy_t, carbohydrate_flow_BufFruits)
    W_Pot_Fruit_1 = fruit_growth(1, last_24_canopy_t) * FRUIT_DEVELOPMENT_STAGES_NUM
    return W_Pot_Fruit_1 * number_flow_BufFruit_1


def carbohydrate_flow_from_buffer_to_fruit_stages(jth: int, carbohydrate_amount_Buf, number_Fruits, canopy_t,
                                                  sum_canopy_t, last_24_canopy_t):
    """
    Equations 9.36
    carbohydrate_flow_BufFruit_j = sum_carbohydrate_flow_BufFruit_conversion_factor * number_fruit_j
                                    * fruit_growth_rate * (carbohydrate_flow_BufFruits - carbohydrate_flow_BufFruit_1)
    Returns: carbohydrates flow from buffer to fruit stage j [mg m^-2 s^-1]
    """
    if jth == 1:
        return carbohydrate_flow_from_buffer_to_first_fruit_stage(carbohydrate_amount_Buf, canopy_t, sum_canopy_t,
                                                                  last_24_canopy_t)

    sum_carbohydrate_flow_BufFruit_conversion_factor = sum_carbohydrate_flow_BufFruit_conversion(number_Fruits, last_24_canopy_t)
    number_fruit_j = number_Fruits[jth]
    fruit_growth_rate_j = fruit_growth(jth, last_24_canopy_t)
    carbohydrate_flow_BufFruits = carbohydrate_flow_from_buffer_to_fruits(carbohydrate_amount_Buf, canopy_t,
                                                                          sum_canopy_t, last_24_canopy_t)
    carbohydrate_flow_BufFruit_1 = carbohydrate_flow_from_buffer_to_fruit_stages(1, carbohydrate_amount_Buf, None,
                                                                                 canopy_t, sum_canopy_t, last_24_canopy_t)
    return sum_carbohydrate_flow_BufFruit_conversion_factor * number_fruit_j \
           * fruit_growth_rate_j * (carbohydrate_flow_BufFruits - carbohydrate_flow_BufFruit_1)


def carbohydrate_flow_from_growth_respiration(carbohydrate_amount_Buf, canopy_t, sum_canopy_t, last_24_canopy_t):
    """
    Equations 9.43
    carbohydrate_flow_BufAir =
    Returns: carbohydrates flow from growth respiration [mg m^-2 s^-1]
    """
    carbohydrate_flow_BufFruits = carbohydrate_flow_from_buffer_to_fruits(carbohydrate_amount_Buf, canopy_t,
                                                                          sum_canopy_t, last_24_canopy_t)
    carbohydrate_flow_BufLeaf = carbohydrate_flow_from_buffer_to_leaves(carbohydrate_amount_Buf, last_24_canopy_t)
    carbohydrate_flow_BufStem = carbohydrate_flow_from_buffer_to_stem(carbohydrate_amount_Buf, last_24_canopy_t)
    return FRUIT_GROWTH_RESPIRATION_COEF * carbohydrate_flow_BufFruits \
           + LEAF_GROWTH_RESPIRATION_COEF * carbohydrate_flow_BufLeaf \
           + STEM_GROWTH_RESPIRATION_COEF * carbohydrate_flow_BufStem


def carbohydrate_flow_from_fruit_maintenance_respiration(carbohydrate_amount_Fruit_j, last_24_canopy_t):
    """
    Equations 9.45
    carbohydrate_flow_FruitAir_j = FRUIT_MAINTENANCE_RESPIRATION_COEF * Q10_M**(0.1*(last_24_canopy_t-25)) * carbohydrate_amount_Fruit_j
                                * (1 - math.exp(-MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF*RELATIVE_GROWTH_RATE))
    Returns: carbohydrates flow from fruit maintenance respiration [mg m^-2 s^-1]
    """
    return FRUIT_MAINTENANCE_RESPIRATION_COEF * Q10_M**(0.1*(last_24_canopy_t-25)) * carbohydrate_amount_Fruit_j \
           * (1 - math.exp(-MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF*RELATIVE_GROWTH_RATE))


def carbohydrate_flow_from_leaf_maintenance_respiration(carbohydrate_amount_Leaf, last_24_canopy_t):
    """
    Equations 9.45
    carbohydrate_flow_LeafAir = LEAF_MAINTENANCE_RESPIRATION_COEF * Q10_M**(0.1*(last_24_canopy_t-25)) * carbohydrate_amount_Leaf \
                                * (1 - math.exp(-MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF*RELATIVE_GROWTH_RATE))
    Returns: carbohydrates flow from leaf maintenance respiration [mg m^-2 s^-1]
    """
    return LEAF_MAINTENANCE_RESPIRATION_COEF * Q10_M**(0.1*(last_24_canopy_t-25)) * carbohydrate_amount_Leaf \
           * (1 - math.exp(-MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF*RELATIVE_GROWTH_RATE))


def carbohydrate_flow_from_stem_maintenance_respiration(carbohydrate_amount_Stem, last_24_canopy_t):
    """
    Equations 9.45
    carbohydrate_flow_StemAir = STEM_MAINTENANCE_RESPIRATION_COEF * Q10_M**(0.1*(last_24_canopy_t-25)) * carbohydrate_amount_Stem \
                                * (1 - math.exp(-MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF*RELATIVE_GROWTH_RATE))
    Returns: carbohydrates flow from stem maintenance respiration [mg m^-2 s^-1]
    """
    return STEM_MAINTENANCE_RESPIRATION_COEF * Q10_M**(0.1*(last_24_canopy_t-25)) * carbohydrate_amount_Stem \
           * (1 - math.exp(-MAINTENANCE_RESPIRATION_FUNCTION_REGRESSION_COEF*RELATIVE_GROWTH_RATE))


def leaf_harvest_rate(carbohydrate_amount_Leaf):
    """
    Equations 9.47, B.5
    carbohydrate_flow_LeafHar = S_ * (carbohydrate_amount_Leaf - maximum_carbohydrates_stored_in_the_leaves)
    Returns: leaf harvest rate [mg m^-2 s^-1]
    """
    max_stored_leaves_carbohydrates = maximum_carbohydrates_stored_in_the_leaves()
    return smoothed_conditional_function(carbohydrate_amount_Leaf, -5e-5, max_stored_leaves_carbohydrates) \
           * (carbohydrate_amount_Leaf - max_stored_leaves_carbohydrates)
