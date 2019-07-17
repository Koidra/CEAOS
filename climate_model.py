import numpy as np
from abc import ABC, abstractmethod
from configs import Config as C
from data_models import Setpoint
from equations import capacities, heatfluxes, lumpedcoverlayers


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
    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        pass


class GreenhouseClimateModel(IndoorClimateModel):

    def __init__(self, greenhouse_config_file):
        super(GreenhouseClimateModel, self).__init__(greenhouse_config_file)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError

    def __canopy_temperature_delta(self, setpoint: Setpoint):
        capCan = capacities.canopy_heat_capacity(C.)


class VerticalFarmClimateModel(IndoorClimateModel):

    def __init__(self, vertical_farm_struct):
        super(VerticalFarmClimateModel, self).__init__(vertical_farm_struct)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError

