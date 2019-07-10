import numpy as np
from abc import ABC, abstractmethod


class IndoorClimateModel(ABC):
    """ An indoor climate model
    Takes a climate control setpoint vector and returns a climate observation vector as output
    """

    # Greenhouse or Vertical farm design constants
    # E.g.: Lamps Par, Light Transmission Ratio, ...
    gh_structure = None

    def __init__(self, greenhouse_struct):
        self.gh_structure = greenhouse_struct

    @abstractmethod
    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        pass


class GreenhouseClimateModel(IndoorClimateModel):

    def __init__(self, greenhouse_struct):
        super(GreenhouseClimateModel, self).__init__(greenhouse_struct)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError


class VerticalFarmClimateModel(IndoorClimateModel):

    def __init__(self, vertical_farm_struct):
        super(VerticalFarmClimateModel, self).__init__(vertical_farm_struct)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError

