import math

from constants import *
from crop.tomato.equations.radiations import PAR_absorbed_by_canopy
from crop.tomato.equations.utils import leaf_area_index


def electron_transport(carbohydrate_amount_Leaf, outdoor_global_rad, canopy_t):
    """
    Equation 9.14
    electron_transport_rate =  (potential_electron_transport_rate + PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR * canopy_PAR_absorbed
                                - sqrt((potential_electron_transport_rate + PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR * canopy_PAR_absorbed)**2
                                        - 4 * ELECTRON_TRANSPORT_RATE_CURVATURE * potential_electron_transport_rate
                                            * PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR * canopy_PAR_absorbed)) \
                                /(2 * ELECTRON_TRANSPORT_RATE_CURVATURE)

    Returns: electron transport rate [µmol {e-} m^-2 s^-1]

    """
    potential_electron_transport_rate = potential_electron_transport(carbohydrate_amount_Leaf, canopy_t)
    canopy_PAR_absorbed = PAR_absorbed_by_canopy(carbohydrate_amount_Leaf, outdoor_global_rad)
    return (potential_electron_transport_rate + PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR * canopy_PAR_absorbed
            - math.sqrt((potential_electron_transport_rate + PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR * canopy_PAR_absorbed) ** 2
                        - 4 * ELECTRON_TRANSPORT_RATE_CURVATURE * potential_electron_transport_rate
                        * PHOTONS_TO_ELECTRONS_CONVERSION_FACTOR * canopy_PAR_absorbed)) \
            / (2 * ELECTRON_TRANSPORT_RATE_CURVATURE)


def potential_electron_transport(carbohydrate_amount_Leaf, canopy_t):
    """
    Equation 9.15
        potential_electron_transport_rate
    = max_electron_transport_rate_at_25
    * math.exp(ACTIVATION_ENERGY_JPOT *
               (reference_canopy_t - REFERENCE_TEMPERATURE_JPOT)/(M_GAS*reference_canopy_t*REFERENCE_TEMPERATURE_JPOT)))
    * (1+math.exp((ENTROPY_TERM_JPOT * REFERENCE_TEMPERATURE_JPOT - DEACTIVATION_ENERGY_JPOT)/(M_GAS*REFERENCE_TEMPERATURE_JPOT)))
    / (1+math.exp((ENTROPY_TERM_JPOT * reference_canopy_t - DEACTIVATION_ENERGY_JPOT)/(M_GAS*reference_canopy_t)))

    Returns: potential electron transport rate[µmol {e-} m^-2 s^-1]
    """
    reference_canopy_t = canopy_t + 273.15
    max_canopy_electron_transport_rate_at_25 = max_canopy_electron_transport_at_25(carbohydrate_amount_Leaf)
    return max_canopy_electron_transport_rate_at_25 \
        * math.exp(ACTIVATION_ENERGY_JPOT
                   * (reference_canopy_t - REFERENCE_TEMPERATURE_JPOT)/(M_GAS*reference_canopy_t*REFERENCE_TEMPERATURE_JPOT)) \
        * (1+math.exp((ENTROPY_TERM_JPOT * REFERENCE_TEMPERATURE_JPOT - DEACTIVATION_ENERGY_JPOT)/(M_GAS*REFERENCE_TEMPERATURE_JPOT))) \
        / (1+math.exp((ENTROPY_TERM_JPOT * reference_canopy_t - DEACTIVATION_ENERGY_JPOT)/(M_GAS*reference_canopy_t)))


def max_canopy_electron_transport_at_25(carbohydrate_amount_Leaf):
    """
    Equation 9.16
    max_electron_transport_rate_at_25 = LAI * MAX_LEAF_ELECTRON_TRANSPORT_RATE
    Returns: maximum rate of electron transport at 25°C for the canopy [µmol {e-} m^-2 s^-1]
    """
    return leaf_area_index(carbohydrate_amount_Leaf) * MAX_LEAF_ELECTRON_TRANSPORT_RATE
