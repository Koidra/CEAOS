from constants import *


def buffer_carbohydrates_amount(crop_inputs: CropInputs):
    """
    Equation 9.1
    carbohydrate_amount_Buf = carbohydrate_mass_flow_AirBuf  − carbohydrate_mass_flow_BufFruit
          − carbohydrate_mass_flow_BufLeaf − carbohydrate_mass_flow_BufStem − carbohydrate_mass_flow_BufAir
    Returns: The evolution of the carbohydrates in the buffer [mg m^-2 s^-1]
    """
    carbohydrate_mass_flow_AirBuf = 0
    carbohydrate_mass_flow_BufFruit = 0
    carbohydrate_mass_flow_BufLeaf = 0
    carbohydrate_mass_flow_BufStem = 0
    carbohydrate_mass_flow_BufAir = 0
    return carbohydrate_mass_flow_AirBuf - carbohydrate_mass_flow_BufFruit \
         - carbohydrate_mass_flow_BufLeaf - carbohydrate_mass_flow_BufStem - carbohydrate_mass_flow_BufAir


def fruit_development_stored_carbohydrates_amount(jth: int, crop_inputs: CropInputs):
    """
    Equation 9.2
    carbohydrate_amount_Fruit_j = carbohydrate_mass_flowBufFruit_j + carbohydrate_mass_flowFruit_jminus_Fruit_j
               − carbohydrate_mass_flow_Fruit_j_jplus - carbohydrate_mass_flow_FruitAir_j
    Returns: Carbohydrates are stored in the fruit development stage j [mg m^-2 s^-1]
    """
    carbohydrate_mass_flowBufFruit_j = 0
    carbohydrate_mass_flowFruit_jminus_Fruit_j = 0
    carbohydrate_mass_flow_Fruit_j_jplus = 0
    carbohydrate_mass_flow_FruitAir_j = 0
    return carbohydrate_mass_flowBufFruit_j + carbohydrate_mass_flowFruit_jminus_Fruit_j \
         - carbohydrate_mass_flow_Fruit_j_jplus - carbohydrate_mass_flow_FruitAir_j


def number_of_fruits(jth: int, crop_inputs: CropInputs):
    """
    Equation 9.3
    number_Fruit{j} = number_flow_Fruit_jminus_Fruit_j - number_flow_Fruit_j_Fruit_jplus
    Returns: The number of fruits in the fruit development stage j  [fruits m^-2 s^-1]
    """
    number_flow_Fruit_jminus_Fruit_j = 0
    number_flow_Fruit_j_Fruit_jplus = 0
    return number_flow_Fruit_jminus_Fruit_j - number_flow_Fruit_j_Fruit_jplus


def leaves_stored_carbohydrates_amount(crop_inputs: CropInputs):
    """
    Equation 9.4
    carbohydrate_amount_Leaf = carbohydrate_mass_BufLeaf - carbohydrate_mass_LeafAir - carbohydrate_mass_LeafHar
    Returns: The carbohydrates stored in the leaves [mg m^-2 s^-1]
    """
    carbohydrate_mass_BufLeaf = 0
    carbohydrate_mass_LeafAir = 0
    carbohydrate_mass_LeafHar = 0
    return carbohydrate_mass_BufLeaf - carbohydrate_mass_LeafAir - carbohydrate_mass_LeafHar


def stem_and_roots_stored_carbohydrates_amount(crop_inputs: CropInputs):
    """
    Equation 9.6
    carbohydrate_amount_Stem = carbohydrate_mass_BufStem - carbohydrate_mass_StemAir
    Returns: The carbohydrates stored in the stem and roots [mg m^-2 s^-1]
    """
    carbohydrate_mass_BufStem = 0
    carbohydrate_mass_StemAir = 0
    return carbohydrate_mass_BufStem - carbohydrate_mass_StemAir


def accumulated_harvested_tomato_dry_matter(crop_inputs: CropInputs):
    """
    Equation 9.7
    dry_matter_Har = CARBOHYDRATE_TO_DRY_MATTER_CONVERSION * carbohydrate_mass_FruitHar
    Returns: THe accumulated harvested tomato dry matter [mg {DM} m^-2 s^-1]
    """
    carbohydrate_mass_FruitHar = 0
    return CARBOHYDRATE_TO_DRY_MATTER_CONVERSION * carbohydrate_mass_FruitHar


def temperature_sum(crop_inputs: CropInputs):
    """
    Equation 9.8
    sum_canopy_t = canopy_t/86400
    Returns: The temperature sum [°C s^-1]
    """
    return crop_inputs.canopy_t/86400


def _24_mean_temperature(crop_inputs: CropInputs):
    """
    Equation 9.9
    _24_canopy_t = 1/DAY_MEAN_TEMP_TIME_CONSTANT * (PROCESS_GAIN * canopy_t - _24_canopy_t)
    Returns: The 24 hour mean canopy temperature [°C s^-1]
    """
    return 1/DAY_MEAN_TEMP_TIME_CONSTANT * (PROCESS_GAIN * crop_inputs.canopy_t - crop_inputs._24_canopy_t)

