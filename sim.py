import numpy as np
from numpy import ndarray as vec
from climate_model import IndoorClimateModel
from crop_model import CropModel


class Simulator(object):
    def __init__(self, climate_model: IndoorClimateModel, crop_model: CropModel):
        self.climate_model = climate_model
        self.crop_model = crop_model
        self._crop_observations = NotImplemented

    def step(self, climate_setpoint: vec, crop_setpoint: vec) -> vec:
        climate_observations = self.climate_model.step(climate_setpoint, self._crop_observations)
        crop_observations = self.crop_model.step(climate_observations, crop_setpoint)
        self._crop_observations = crop_observations
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def render(self, mode='human'):
        pass
