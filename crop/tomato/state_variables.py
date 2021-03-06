from ...data_models import ClimateStates
from .crop_model import CropStates
from .carbohydrate_flows import *
from .fruit_flow import fruit_flow_through_fruit_development_stage


def buffer_carbohydrates_amount(crop_states: CropStates, climate_states: ClimateStates):
    """
    Equation 9.1
    carbohydrate_amount_Buf = carbohydrate_flow_AirBuf − carbohydrate_flow_BufFruit
          − carbohydrate_flow_BufLeaf − carbohydrate_flow_BufStem − carbohydrate_flow_BufAir
    Returns: The evolution of the carbohydrates in the buffer [mg m^-2 s^-1]
    """
    carbohydrate_flow_AirBuf = net_photosynthesis_rate(crop_states.carbohydrate_amount_Buf,
                                                       crop_states.carbohydrate_amount_Leaf, climate_states.co2_Air,
                                                       climate_states.t_Canopy,
                                                       climate_states.PAR_Canopy)
    carbohydrate_flow_BufFruit = carbohydrate_flow_from_buffer_to_fruits(crop_states.carbohydrate_amount_Buf,
                                                                         climate_states.t_Canopy,
                                                                         crop_states.sum_canopy_t,
                                                                         crop_states.last_24_canopy_t)
    carbohydrate_flow_BufLeaf = carbohydrate_flow_from_buffer_to_leaves(crop_states.carbohydrate_amount_Buf,
                                                                        crop_states.last_24_canopy_t)
    carbohydrate_flow_BufStem = carbohydrate_flow_from_buffer_to_stem(crop_states.carbohydrate_amount_Buf,
                                                                      crop_states.last_24_canopy_t)
    carbohydrate_flow_BufAir = carbohydrate_flow_from_growth_respiration(crop_states.carbohydrate_amount_Buf,
                                                                         climate_states.t_Canopy,
                                                                         crop_states.sum_canopy_t,
                                                                         crop_states.last_24_canopy_t)
    return carbohydrate_flow_AirBuf - carbohydrate_flow_BufFruit \
           - carbohydrate_flow_BufLeaf - carbohydrate_flow_BufStem - carbohydrate_flow_BufAir


def fruit_development_stored_carbohydrates_amount(jth: int, crop_inputs: CropStates, climate_inputs: ClimateStates):
    """
    Equation 9.2
    carbohydrate_amount_Fruit_j = carbohydrate_flow_BufFruit_j + carbohydrate_flow_Fruit_jminus_Fruit_j
               − carbohydrate_flow_Fruit_j_jplus - carbohydrate_flow_FruitAir_j
    Returns: Carbohydrates are stored in the fruit development stage j [mg m^-2 s^-1]
    """
    carbohydrate_flow_BufFruit_j = carbohydrate_flow_from_buffer_to_fruit_stages(jth,
                                                                                 crop_inputs.carbohydrate_amount_Buf,
                                                                                 crop_inputs.number_Fruits,
                                                                                 climate_inputs.t_Canopy,
                                                                                 crop_inputs.sum_canopy_t,
                                                                                 crop_inputs.last_24_canopy_t)
    carbohydrate_flow_Fruit_jminus_Fruit_j = carbohydrate_flow_through_fruit_stages(jth - 1,
                                                                                    crop_inputs.carbohydrate_amount_Fruits,
                                                                                    crop_inputs.last_24_canopy_t)
    carbohydrate_flow_Fruit_j_Fruit_jplus = carbohydrate_flow_through_fruit_stages(jth,
                                                                                   crop_inputs.carbohydrate_amount_Fruits,
                                                                                   crop_inputs.last_24_canopy_t)
    carbohydrate_flow_FruitAir_j = carbohydrate_flow_from_fruit_maintenance_respiration(
        crop_inputs.carbohydrate_amount_Fruits[jth],
        crop_inputs.last_24_canopy_t)
    return carbohydrate_flow_BufFruit_j + carbohydrate_flow_Fruit_jminus_Fruit_j \
           - carbohydrate_flow_Fruit_j_Fruit_jplus - carbohydrate_flow_FruitAir_j


def number_of_fruits(jth: int, crop_states: CropStates):
    """
    Equation 9.3
    number_Fruit_j = number_flow_Fruit_jminus_Fruit_j - number_flow_Fruit_j_Fruit_jplus
    Returns: The number of fruits in the fruit development stage j  [fruits m^-2 s^-1]

    Args:
        crop_states:
    """
    number_flow_Fruit_jminus_Fruit_j = fruit_flow_through_fruit_development_stage(jth - 1, crop_states.number_Fruits,
                                                                                  crop_states.sum_canopy_t,
                                                                                  crop_states.last_24_canopy_t)
    number_flow_Fruit_j_Fruit_jplus = fruit_flow_through_fruit_development_stage(jth, crop_states.number_Fruits,
                                                                                 crop_states.sum_canopy_t,
                                                                                 crop_states.last_24_canopy_t)
    return number_flow_Fruit_jminus_Fruit_j - number_flow_Fruit_j_Fruit_jplus


def leaves_stored_carbohydrates_amount(crop_states: CropStates):
    """
    Equation 9.4
    carbohydrate_amount_Leaf = carbohydrate_flow_BufLeaf - carbohydrate_flow_LeafAir - carbohydrate_flow_LeafHar
    Returns: The carbohydrates stored in the leaves [mg m^-2 s^-1]
    """
    carbohydrate_flow_BufLeaf = carbohydrate_flow_from_buffer_to_leaves(crop_states.carbohydrate_amount_Buf,
                                                                        crop_states.last_24_canopy_t)
    carbohydrate_flow_LeafAir = carbohydrate_flow_from_leaf_maintenance_respiration(crop_states.carbohydrate_amount_Buf,
                                                                                    crop_states.last_24_canopy_t)
    carbohydrate_flow_LeafHar = leaf_harvest_rate(crop_states.carbohydrate_amount_Leaf)
    return carbohydrate_flow_BufLeaf - carbohydrate_flow_LeafAir - carbohydrate_flow_LeafHar


def stem_and_roots_stored_carbohydrates_amount(crop_states: CropStates):
    """
    Equation 9.6
    carbohydrate_amount_Stem = carbohydrate_flow_BufStem - carbohydrate_flow_StemAir
    Returns: The carbohydrates stored in the stem and roots [mg m^-2 s^-1]
    """
    carbohydrate_flow_BufStem = carbohydrate_flow_from_buffer_to_stem(crop_states.carbohydrate_amount_Buf,
                                                                      crop_states.last_24_canopy_t)
    carbohydrate_flow_StemAir = carbohydrate_flow_from_stem_maintenance_respiration(crop_states.carbohydrate_amount_Buf,
                                                                                    crop_states.last_24_canopy_t)
    return carbohydrate_flow_BufStem - carbohydrate_flow_StemAir


def accumulated_harvested_tomato_dry_matter(crop_states: CropStates):
    """
    Equation 9.7
    dry_matter_Har = CARBOHYDRATE_TO_DRY_MATTER_CONVERSION * carbohydrate_flow_FruitHar
    Returns: THe accumulated harvested tomato dry matter [mg {DM} m^-2 s^-1]
    """
    carbohydrate_flow_FruitHar = carbohydrate_flow_through_fruit_stages(FRUIT_DEVELOPMENT_STAGES_NUM,
                                                                        crop_states.carbohydrate_amount_Fruits,
                                                                        crop_states.last_24_canopy_t)
    return CARBOHYDRATE_TO_DRY_MATTER_CONVERSION * carbohydrate_flow_FruitHar


def temperature_sum(climate_states: ClimateStates):
    """
    Equation 9.8
    sum_canopy_t = t_Canopy/86400
    Returns: The temperature sum [°C s^-1]
    """
    return climate_states.t_Canopy / 86400


def _24_mean_temperature(crop_states: CropStates, climate_states: ClimateStates):
    """
    Equation 9.9
    last_24_canopy_t = 1/DAY_MEAN_TEMP_TIME_CONSTANT * (PROCESS_GAIN * t_Canopy - last_24_canopy_t)
    Returns: The 24 hour mean canopy temperature [°C s^-1]
    """
    return 1 / DAY_MEAN_TEMP_TIME_CONSTANT * (PROCESS_GAIN * climate_states.t_Canopy - crop_states.last_24_canopy_t)
