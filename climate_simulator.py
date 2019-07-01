from abc import ABC, abstractmethod

class ClimateSimulator(ABC):
    # a climate simulator takes a climate control setpoint vector and returns a climate obersation vector as output
    def step(setpoint: np.ndarray)

class GreenhouseClimateSimulator(ClimateSimulator):
    def step(setpoint: np.ndarray):
        raise NotImplementedError


class VerticalFarmClimateSimulator(ClimateSimulator):
    def step(setpoint: np.ndarray):
        raise NotImplementedError
