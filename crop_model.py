from abc import ABC, abstractmethod

import numpy as np


class CropModel(ABC):
    @abstractmethod
    def step(self, climate_obs: np.ndarray, setpoint: np.ndarray):
        pass
