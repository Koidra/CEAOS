# CEA Open Simulator (CEAOS)
CEAOS is a dynamic simulator for greenhouses and vertical farms. With this simulator, we hope to accelerate the CEA (Controlled Environment Agriculture) research.

## Concepts  
### Climate Model  
This model estimates the indoor climate based on the outdoor climate, the greenhouse design and the growth of crop (from crop model). For example, these indoor attributes will be estimated by this model: temperature, vapour pressure, CO2,... The estimated indoor climate will be used as the input for the crop model.  

### Crop Model  
This model describes the effects of indoor climate and greenhouse design on the crop. In other words, this model describes crop as a function of indoor climate.  

Both the climate model and crop model should consist of a set of differential equations and to allow the use of the ordinary differential equation solver.


## Interface  
The simulator inherits OpenAI Gym [Env](https://github.com/openai/gym/blob/master/gym/core.py) interface. The following is two main methods of the interface:  
 - `reset(self)`: Reset the simulator's state. Returns initial observations.
 - `step(self, setpoint)`: Step the simulator by one step, return observations.  
 The observations are sensor values, indoor climate values, etc, at the next time step.

The detail to interact with the simulator at the code level is by  
1. First instantiating an environment object with certain settings, this describes the greenhouse design or vertical farm design. E.g.
```python
sim = Greenhouse(crop_name, greenhouse_structural_properties, region, begin_time, end_time, time_step=5, ...)
```

2. Then at each step, call the `step` function to forward the simulator one time step (e.g. 5 minutes) with the setpoint provided by the controller
```python
observations = sim.step(setpoints)
```
The result of `step` is the observations, e.g. sensor values, at the next time step.

## Modular Software Architecture
![CEAOS Architecture](https://i.imgur.com/v1IcJXg.png)

When the `step` method of simulator is called, it references to both `climate model` and `crop model` to estimate coefficients and solve differential equations to get the observations at the current time step. Climate model and crop model share each other information to complete the task. 

Go into details, there are two climate models (Greenhouse model, vertical farm model) and many types of crop model (E.g. Lettuce crop, cucumber crop or tomato crop). So it's required the object-oriented software design, also each inherited class of `crop model` and `climate model` must use the same interface.