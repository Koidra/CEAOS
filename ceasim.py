import numpy as np
from gym import Env
from climate_model import IndoorClimateModel
from crop_model import CropModel


class CEASim(Env):

    def __init__(self, climate_model: IndoorClimateModel, crop_model: CropModel):
        self.climate_model = climate_model
        self.crop_model = crop_model
        self._crop_obs = NotImplemented

    def step_(self, climate_setpoint: np.ndarray, crop_setpoint: np.ndarray) -> np.ndarray:
        climate_obs = self.climate_model.step(climate_setpoint, self._crop_obs)
        crop_obs = self.crop_model.step(climate_obs, crop_setpoint)

        self._crop_obs = crop_obs
        raise NotImplementedError

    # Implement gym's interface
    def step(self, setpoint: np.ndarray):
        # TODO: Call self._step
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def render(self, mode='human'):
        pass
