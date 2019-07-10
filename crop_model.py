import numpy as np
from abc import ABC, abstractmethod


class CropModel(ABC):
    @abstractmethod
    def step(self, climate_obs: np.ndarray, setpoint: np.ndarray):
        pass
