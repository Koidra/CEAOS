from abc import ABC, abstractmethod

import numpy as np

from data_models import Setpoints


class IndoorClimateModel(ABC):
    """ An indoor climate model
    Takes a climate control setpoint vector and returns a climate observation vector as output
    """

    # Greenhouse or Vertical farm design constants
    # E.g.: Lamps Par, Light Transmission Ratio, ...
    greenhouse_design = None

    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def step(self, crop_observations: np.ndarray, setpoint: Setpoints):
        pass


class GreenhouseClimateModel(IndoorClimateModel):

    def __init__(self, greenhouse_config_file):
        super(GreenhouseClimateModel, self).__init__(greenhouse_config_file)

    def step(self, crop_observations: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError


class VerticalFarmClimateModel(IndoorClimateModel):

    def __init__(self, vertical_farm_struct):
        super(VerticalFarmClimateModel, self).__init__(vertical_farm_struct)

    def step(self, crop_observations: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError
