"""The implementation of Lumped Cover Layers' equations

Based on section 8.4
r_: The reflection coefficient of the layer
t_: The transmission coefficient of the layer

The default model contains three cover layers, i.e.
    - A movable outdoor shading screen (ShScr)
    - The greenhouse roof (Rf)
    - A movable indoor thermal screen (ThScr).
"""
import math

from climate.heat_fluxes import canopy_virtual_NIR_reflection_coefficient, \
    lumped_cover_virtual_NIR_transmission_coefficients, canopy_virtual_NIR_transmission_coefficient, \
    floor_virtual_NIR_transmission_coefficients
from coefficients import Coefficients
from data_models import Setpoints, ClimateStates


def double_layer_cover_transmission_coefficient(transmission_coef_1, transmission_coef_2, reflection_coef_1_Dn, reflection_coef_2_Up) -> float:
    """
    The transmission coefficient of a double layer cover
    Equation 8.14, A4 [2]
    :param float transmission_coef_1: the transmission coefficients of the first layer
    :param float transmission_coef_2: the transmission coefficients of the second layer
    :param float reflection_coef_1_Dn: the reflection coefficients of the first layer
    :param float reflection_coef_2_Up: the reflection coefficients of the second layer
    :return: The transmission coefficient
    """
    return (transmission_coef_1 * transmission_coef_2) / (1 - reflection_coef_1_Dn * reflection_coef_2_Up)


def double_layer_cover_up_reflection_coefficient(transmission_coef_1, reflection_coef_1_Up,
                                                 reflection_coef_1_Down, reflection_coef_2_Up) -> float:
    """
    The reflection coefficient of a double layer cover
    Equation 8.15
    :param float transmission_coef_1: the transmission coefficients of the first layer
    :param float reflection_coef_1_Up: the reflection coefficients of the first layer
    :param float reflection_coef_1_Down: the reflection coefficients of the first layer
    :param float reflection_coef_2_Up: the reflection coefficients of the second layer
    :return: The reflection coefficient
    """

    return reflection_coef_1_Up + ((transmission_coef_1**2) * reflection_coef_2_Up) \
                                    / (1 - reflection_coef_1_Down * reflection_coef_2_Up)


def double_layer_cover_down_reflection_coefficient(transmission_coef_2, reflection_coef_1_Down,
                                                   reflection_coef_2_Up, reflection_coef_2_Down) -> float:
    """
    The reflection coefficient of a double layer cover
    Equation 8.15
    :param float transmission_coef_2: the transmission coefficients of the second layer
    :param float reflection_coef_1_Down: the reflection coefficients of the first layer
    :param float reflection_coef_2_Up: the reflection coefficients of the second layer
    :param float reflection_coef_2_Down: the reflection coefficients of the second layer
    :return: The reflection coefficient
    """

    return reflection_coef_2_Down + ((transmission_coef_2**2) * reflection_coef_1_Down) \
                                    / (1 - reflection_coef_1_Down * reflection_coef_2_Up)


def absorption_coefficient(transmission_coef, reflection_coef):
    """The absorption coefficient

    :param transmission_coef: the transmission coefficient
    :param reflection_coef: the reflection coefficient
    :return: The absorption coefficient
    """
    return 1 - (transmission_coef + reflection_coef)


def lumped_cover_heat_capacity():
    # Equation 8.18
    mean_greenhouse_cover_slope = Coefficients.Construction.mean_greenhouse_cover_slope
    roof_thickness = Coefficients.Roof.roof_thickness
    roof_density = Coefficients.Roof.roof_density
    c_p_Rf = Coefficients.Roof.c_p_Rf
    return math.cos(mean_greenhouse_cover_slope)*(roof_thickness*roof_density*c_p_Rf)


def lumped_cover_conductive_heat_flux():
    # Equation 8.19
    roof_thickness = Coefficients.Roof.roof_thickness
    roof_heat_conductivity = Coefficients.Roof.roof_heat_conductivity
    return (roof_thickness/roof_heat_conductivity)**-1


def lumped_cover_NIR_coefficients(setpoints: Setpoints):
    shScr_NIR_transmission_coef = 1 - setpoints.U_ShScr \
                                      * (1 - Coefficients.Shadowscreen.shScr_NIR_transmission_coefficient)  # line 155 / setGlParams / GreenLight
    shScr_NIR_reflection_coef = setpoints.U_ShScr * Coefficients.Shadowscreen.shScr_NIR_reflection_coefficient  # line 152 / setGlParams / GreenLight
    roof_NIR_transmission_coef = setpoints.U_Roof * Coefficients.Roof.roof_NIR_transmission_coefficient
    roof_NIR_reflection_coef = setpoints.U_Roof * Coefficients.Roof.roof_NIR_reflection_coefficient
    thScr_NIR_transmission_coef = 1 - setpoints.U_ThScr \
                                  * (1 - Coefficients.Thermalscreen.thScr_NIR_transmission_coefficient)
    thScr_NIR_reflection_coef = setpoints.U_ThScr * Coefficients.Thermalscreen.thScr_NIR_reflection_coefficient
    blScr_NIR_transmission_coef = 1 - setpoints.U_BlScr * (1 - Coefficients.Blackoutscreen.blScr_NIR_transmission_coef)
    blScr_NIR_reflection_coef = setpoints.U_BlScr * Coefficients.Blackoutscreen.blScr_NIR_reflection_coef
    lamp_NIR_transmission_coef = setpoints.U_Roof * Coefficients.Lamp.lamp_NIR_transmission_coef
    lamp_NIR_reflection_coef = setpoints.U_Roof * Coefficients.Lamp.lamp_NIR_reflection_coef

    # NIR reflection coefficient of the roof and thermal screen
    roof_thScr_NIR_transmission_coef = double_layer_cover_transmission_coefficient(roof_NIR_transmission_coef,
                                                                                   thScr_NIR_transmission_coef,
                                                                                   roof_NIR_reflection_coef,
                                                                                   thScr_NIR_reflection_coef)
    roof_thScr_NIR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(roof_NIR_transmission_coef,
                                                                                     roof_NIR_reflection_coef,
                                                                                     roof_NIR_reflection_coef,
                                                                                     thScr_NIR_reflection_coef)
    roof_thScr_NIR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(thScr_NIR_transmission_coef,
                                                                                         roof_NIR_reflection_coef,
                                                                                         thScr_NIR_reflection_coef,
                                                                                         thScr_NIR_reflection_coef)

    # Vanthoor NIR reflection coefficient of the lumped cover
    vanilla_cover_NIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_NIR_transmission_coef,
                                                                                      roof_thScr_NIR_transmission_coef,
                                                                                      shScr_NIR_reflection_coef,
                                                                                      roof_thScr_NIR_reflection_up_coef)
    vanilla_cover_NIR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(shScr_NIR_transmission_coef,
                                                                                        shScr_NIR_reflection_coef,
                                                                                        shScr_NIR_reflection_coef,
                                                                                        roof_thScr_NIR_reflection_up_coef)
    vanilla_cover_NIR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(
        roof_thScr_NIR_transmission_coef,
        shScr_NIR_reflection_coef,
        roof_thScr_NIR_reflection_up_coef,
        roof_thScr_NIR_reflection_down_coef)

    # NIR reflection coefficient of the lumped cover with new blackout screen
    blScr_cover_NIR_transmission_coef = double_layer_cover_transmission_coefficient(vanilla_cover_NIR_transmission_coef,
                                                                                    blScr_NIR_transmission_coef,
                                                                                    vanilla_cover_NIR_reflection_down_coef,
                                                                                    blScr_NIR_reflection_coef)
    blScr_cover_NIR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(
        vanilla_cover_NIR_transmission_coef,
        vanilla_cover_NIR_reflection_up_coef,
        vanilla_cover_NIR_reflection_down_coef,
        blScr_NIR_reflection_coef)
    blScr_cover_NIR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(blScr_NIR_transmission_coef,
                                                                                          vanilla_cover_NIR_reflection_down_coef,
                                                                                          blScr_NIR_reflection_coef,
                                                                                          blScr_NIR_reflection_coef)

    # NIR reflection coefficient of the lumped cover with new blackout screen and lamp
    cover_NIR_transmission_coef = double_layer_cover_transmission_coefficient(blScr_cover_NIR_transmission_coef,
                                                                              lamp_NIR_transmission_coef,
                                                                              blScr_cover_NIR_reflection_down_coef,
                                                                              lamp_NIR_reflection_coef)
    cover_NIR_reflection_coef = double_layer_cover_up_reflection_coefficient(blScr_cover_NIR_transmission_coef,
                                                                             blScr_cover_NIR_reflection_up_coef,
                                                                             blScr_cover_NIR_reflection_down_coef,
                                                                             lamp_NIR_reflection_coef)
    return cover_NIR_transmission_coef, cover_NIR_reflection_coef


def lumped_cover_PAR_coefficients(setpoints: Setpoints):
    shScr_PAR_transmission_coef = 1 - setpoints.U_ShScr \
                                      * (1 - Coefficients.Shadowscreen.shScr_PAR_transmission_coefficient)  # line 155 / setGlParams / GreenLight
    shScr_PAR_reflection_coef = setpoints.U_ShScr * Coefficients.Shadowscreen.shScr_PAR_reflection_coefficient  # line 152 / setGlParams / GreenLight
    roof_PAR_transmission_coef = setpoints.U_Roof * Coefficients.Roof.roof_PAR_transmission_coefficient
    roof_PAR_reflection_coef = setpoints.U_Roof * Coefficients.Roof.roof_PAR_reflection_coefficient
    thScr_PAR_transmission_coef = 1 - setpoints.U_ThScr * (
                1 - Coefficients.Thermalscreen.thScr_PAR_transmission_coefficient)
    thScr_PAR_reflection_coef = setpoints.U_ThScr * Coefficients.Thermalscreen.thScr_PAR_reflection_coefficient
    blScr_PAR_transmission_coef = 1 - setpoints.U_BlScr * (1 - Coefficients.Blackoutscreen.blScr_PAR_transmission_coef)
    blScr_PAR_reflection_coef = setpoints.U_BlScr * Coefficients.Blackoutscreen.blScr_PAR_reflection_coef
    lamp_PAR_transmission_coef = setpoints.U_Roof * Coefficients.Lamp.lamp_PAR_transmission_coef
    lamp_PAR_reflection_coef = setpoints.U_Roof * Coefficients.Lamp.lamp_PAR_reflection_coef
    # PAR reflection coefficient of the roof and thermal screen
    roof_thScr_PAR_transmission_coef = double_layer_cover_transmission_coefficient(roof_PAR_transmission_coef,
                                                                                   thScr_PAR_transmission_coef,
                                                                                   roof_PAR_reflection_coef,
                                                                                   thScr_PAR_reflection_coef)
    roof_thScr_PAR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(roof_PAR_transmission_coef,
                                                                                     roof_PAR_reflection_coef,
                                                                                     roof_PAR_reflection_coef,
                                                                                     thScr_PAR_reflection_coef)
    roof_thScr_PAR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(thScr_PAR_transmission_coef,
                                                                                         roof_PAR_reflection_coef,
                                                                                         thScr_PAR_reflection_coef,
                                                                                         thScr_PAR_reflection_coef)

    # Vanthoor PAR reflection coefficient of the lumped cover
    vanilla_cover_PAR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_PAR_transmission_coef,
                                                                                      roof_thScr_PAR_transmission_coef,
                                                                                      shScr_PAR_reflection_coef,
                                                                                      roof_thScr_PAR_reflection_up_coef)
    vanilla_cover_PAR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(shScr_PAR_transmission_coef,
                                                                                        shScr_PAR_reflection_coef,
                                                                                        shScr_PAR_reflection_coef,
                                                                                        roof_thScr_PAR_reflection_up_coef)
    vanilla_cover_PAR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(
        roof_thScr_PAR_transmission_coef,
        shScr_PAR_reflection_coef,
        roof_thScr_PAR_reflection_up_coef,
        roof_thScr_PAR_reflection_down_coef)

    # PAR reflection coefficient of the lumped cover with new blackout screen
    blScr_cover_PAR_transmission_coef = double_layer_cover_transmission_coefficient(vanilla_cover_PAR_transmission_coef,
                                                                                    blScr_PAR_transmission_coef,
                                                                                    vanilla_cover_PAR_reflection_down_coef,
                                                                                    blScr_PAR_reflection_coef)
    blScr_cover_PAR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(
        vanilla_cover_PAR_transmission_coef,
        vanilla_cover_PAR_reflection_up_coef,
        vanilla_cover_PAR_reflection_down_coef,
        blScr_PAR_reflection_coef)
    blScr_cover_PAR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(blScr_PAR_transmission_coef,
                                                                                          vanilla_cover_PAR_reflection_down_coef,
                                                                                          blScr_PAR_reflection_coef,
                                                                                          blScr_PAR_reflection_coef)

    # PAR reflection coefficient of the lumped cover with new blackout screen and lamp
    cover_PAR_transmission_coef = double_layer_cover_transmission_coefficient(blScr_cover_PAR_transmission_coef,
                                                                              lamp_PAR_transmission_coef,
                                                                              blScr_cover_PAR_reflection_down_coef,
                                                                              lamp_PAR_reflection_coef)
    cover_PAR_reflection_coef = double_layer_cover_up_reflection_coefficient(blScr_cover_PAR_transmission_coef,
                                                                             blScr_cover_PAR_reflection_up_coef,
                                                                             blScr_cover_PAR_reflection_down_coef,
                                                                             lamp_PAR_reflection_coef)
    return cover_PAR_transmission_coef, cover_PAR_reflection_coef


def cover_canopy_floor_NIR_coefficients(states: ClimateStates, setpoints: Setpoints):
    cover_NIR_transmission_coef, cover_NIR_reflection_coef = lumped_cover_NIR_coefficients(setpoints)

    virtual_NIR_reflection_canopy_coef = canopy_virtual_NIR_reflection_coefficient(states)
    floor_NIR_reflection_coef = Coefficients.Floor.floor_NIR_reflection_coefficient

    virtual_cover_NIR_transmission_coef = lumped_cover_virtual_NIR_transmission_coefficients(cover_NIR_reflection_coef)
    virtual_canopy_NIR_transmission_coef = canopy_virtual_NIR_transmission_coefficient(states)
    virtual_floor_NIR_transmission_coef = floor_virtual_NIR_transmission_coefficients()

    # NIR transmission coefficient of the cover and canopy
    cover_canopy_NIR_transmission_coef = double_layer_cover_transmission_coefficient(
        virtual_cover_NIR_transmission_coef,
        virtual_canopy_NIR_transmission_coef,
        cover_NIR_reflection_coef,
        virtual_NIR_reflection_canopy_coef)  # line 380 / setGlAux / GreenLight

    # NIR reflection coefficient of the cover and canopy
    cover_canopy_NIR_reflection_up_coef = double_layer_cover_up_reflection_coefficient(
        virtual_cover_NIR_transmission_coef,
        cover_NIR_reflection_coef,
        cover_NIR_reflection_coef,
        virtual_NIR_reflection_canopy_coef)  # line 383, 386 / setGlAux / GreenLight

    cover_canopy_NIR_reflection_down_coef = double_layer_cover_down_reflection_coefficient(
        virtual_canopy_NIR_transmission_coef,
        cover_NIR_reflection_coef,
        virtual_NIR_reflection_canopy_coef,
        virtual_NIR_reflection_canopy_coef)

    # NIR transmission coefficient of the cover, canopy and floor
    cover_canopy_floor_NIR_transmission_coef = double_layer_cover_transmission_coefficient(
        cover_canopy_NIR_transmission_coef,
        virtual_floor_NIR_transmission_coef,
        cover_canopy_NIR_reflection_down_coef,
        floor_NIR_reflection_coef)  # line 389 / setGlAux / GreenLight

    cover_canopy_floor_NIR_reflection_coef = double_layer_cover_up_reflection_coefficient(
        cover_canopy_NIR_transmission_coef,
        cover_canopy_NIR_reflection_up_coef,
        cover_canopy_NIR_reflection_down_coef,
        floor_NIR_reflection_coef)

    return cover_canopy_floor_NIR_transmission_coef, cover_canopy_floor_NIR_reflection_coef


def cover_FIR_coefficients(setpoints: Setpoints):
    shScr_FIR_transmission_coef = 1 - setpoints.U_ShScr*(1 - Coefficients.Shadowscreen.shScr_FIR_transmission_coefficient)
    shScr_FIR_reflection_coef = Coefficients.Shadowscreen.shScr_FIR_reflection_coefficient
    roof_FIR_transmission_coef = Coefficients.Roof.roof_FIR_transmission_coefficient
    roof_FIR_reflection_coef = Coefficients.Roof.roof_FIR_reflection_coefficient
    cover_FIR_transmission_coef = double_layer_cover_transmission_coefficient(shScr_FIR_transmission_coef,
                                                                              roof_FIR_transmission_coef,
                                                                              shScr_FIR_reflection_coef,
                                                                              roof_FIR_reflection_coef)  # line 255 / setGlAux / GreenLight
    cover_FIR_reflection_coef = double_layer_cover_up_reflection_coefficient(
        shScr_FIR_transmission_coef,
        shScr_FIR_reflection_coef,
        shScr_FIR_reflection_coef,
        roof_FIR_reflection_coef)  # line 260 / setGlAux / GreenLight
    return cover_FIR_transmission_coef, cover_FIR_reflection_coef
