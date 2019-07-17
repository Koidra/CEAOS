"""The implementation of Lumped Cover Layers' equations

Based on section 8.4
r_: The reflection coefficient of the layer
t_: The transmission coefficient of the layer

The default model contains four cover layers, i.e.
a movable outdoor shading screen (ShScr),
a semi-permanent shading screen (ShScrPer),
the greenhouse roof (Rf) and
a movable indoor thermal screen (ThScr).
"""


def transmission_coefficient(t1, t2, r1, r2) -> float:
    """
    The transmission coefficient
    Equation 8.14
    :param float t1: the transmission coefficients of the first layer
    :param float t2: the transmission coefficients of the second layer
    :param float r1: the reflection coefficients of the first layer
    :param float r2: the reflection coefficients of the second layer
    :return: The transmission coefficient
    """

    return (t1 * t2) / (1 - r1 * r2)


def reflection_coefficient(t1, r1, r2) -> float:
    """
    The reflection coefficient
    Equation 8.15
    :param float t1: the transmission coefficients of the first layer
    :param float r1: the reflection coefficients of the first layer
    :param float r2: the reflection coefficients of the second layer
    :return: The reflection coefficient
    """

    return r1 + (t1 * t1 * r2) / (1 - r1 * r2)


def transmission_coefficient_with_control(u1, u2, t1, t2, r1, r2) -> float:
    """
    The transmission coefficient with control
    Equation 8.16
    :param float u1: The control of layer 1
    :param float u2: The control of layer 2
    :param float t1: the transmission coefficients of the first layer
    :param float t2: the transmission coefficients of the second layer
    :param float r1: the reflection coefficients of the first layer
    :param float r2: the reflection coefficients of the second layer
    :return: The transmission coefficient with control
    """

    t1_ = (1 - u1 * (1 - t1))
    t2_ = 1 - u2 * (1 - t2)
    r1_ = u1 * r1
    r2_ = u2 * r2

    return transmission_coefficient(t1, t2, r1, r2)


def reflection_coefficient_with_control(u1, u2, t1, t2, r1, r2) -> float:
    """
    The reflection coefficient with control
    Equation 8.17
    :param float u1: The control of layer 1
    :param float u2: The control of layer 2
    :param float t1: the transmission coefficients of the first layer
    :param float t2: the transmission coefficients of the second layer
    :param float r1: the reflection coefficients of the first layer
    :param float r2: the reflection coefficients of the second layer
    :return: The reflection coefficient with control
    """
    t1_ = (1 - u1 * (1 - t1))
    t2_ = 1 - u2 * (1 - t2)
    r1_ = u1 * r1
    r2_ = u2 * r2

    return reflection_coefficient(t1, r1, r2)


def absorption_coefficient(t, r):
    """The absorption coefficient

    :param t: the transmission coefficient
    :param r: the reflection coefficient
    :return: The absorption coefficient
    """
    return 1 - (t + r)
