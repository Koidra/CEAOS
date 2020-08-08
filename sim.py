import numpy as np
import np.ndarray as vec
from climate_model import IndoorClimateModel
from crop_model import CropModel


class Simulator(object):
    def __init__(self, climate_model: IndoorClimateModel, crop_model: CropModel):
        self.climate_model = climate_model
        self.crop_model = crop_model
        self._crop_obs = NotImplemented

    def step(self, climate_setpoint: vec, crop_setpoint: vec) -> vec:
        climate_obs = self.climate_model.step(climate_setpoint, self._crop_obs)
        crop_obs = self.crop_model.step(climate_obs, crop_setpoint)

        self._crop_obs = crop_obs
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def render(self, mode='human'):
        pass
