from abc import ABC, abstractmethod
import numpy as np
from gym import Env
import climate_simulator, crop_simulator


class CEASim(Env):
    def __init__(structure, crop):
        self.climatesim = NotImplemented
        self.cropsim = NotImplemented
        self._state = NotImplemented

    def step_(self, climate_setpoint: np.ndarray, crop_setpoint: np.ndarray) -> np.ndarray:
        climate_obs = self.climatesim.step(climate_setpoint, self._state)
        crop_obs = self.climatesim.step(...)
        self._state = #update state

        raise NotImplementedError()


    # Implement gym's interface
    def step(self, setpoint: np.ndarray):
        return self.step_(setpoint split)
