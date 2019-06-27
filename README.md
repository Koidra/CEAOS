# CEA Open Simulator (CEAOS)
CEAOS is a dynamic simulator for greenhouses and vertical farms. With this simulator, we hope to accelerate the CEA (Controlled Environment Agriculture) research.

## Interface
The way to interact with the simulator at the code level is by

1. First instantiating an environment object with certain settings, e.g.
```python
sim = Greenhouse(crop_name, greenhouse_structural_properties, region, begin_time, end_time, time_step=5, ...)
```

2. Then at each step, call the `step` function to forward the simulator one time step (e.g. 5 minutes) with the setpoint provided by the controller
```python
observations = sim.step(setpoints)
```
The result of `step` is the observations, e.g. sensor values, at the next time step.

## Modular Software Architecture
![CEAOS Architecture](https://i.snag.gy/UV28I4.jpg)
