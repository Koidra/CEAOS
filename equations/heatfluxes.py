"""The implementation of Heat fluxes' equations

Based on section 8.6
"""
import numpy as np


# 8.6.1 Global, PAR and NIR heat fluxes


def R_PARGh(h_GlobAir, t_CovPar, h_GlobPAR, i_Glob):
    """The PAR above the canopy
    Equation 8.28

    :param float h_GlobAir: the ratio of the global radiation
    :param float t_CovPar: the PAR transmission coefficient of the greenhouse cover
    :param float h_GlobPAR: the ratio between PAR and the global radiation
    :param float i_Glob: the outside global radiation
    :return: the PAR above the canopy [W*(m^-2)]
    """
    return (1 - h_GlobAir) * t_CovPar * h_GlobPAR * i_Glob


def R_PAR_SunCan(R_PARGh_, r_CanPAR, r_FlrPAR, K1_PAR, K2_PAR, LAI):
    """The PAR absorbed by the canopy
    Equation 8.26

    :param float R_PARGh_: The PAR above the canopy, from function `R_PARGh`
    :param float r_CanPAR: the reflection coefficient of the canopy for PAR
    :param float r_FlrPAR: the reflection coefficient of the greenhouse floor for PAR
    :param float K1_PAR: the extinction coefficient of the canopy for PAR
    :param float K2_PAR: the extinction coefficient for PAR that is reflected from the floor to the canopy
    :param float LAI: Leaf Area Index
    :return: The PAR absorbed by the canopy [W*(m^-2)]
    """
    R_PAR_SunCan_dn = R_PARGh_ * (1 - r_CanPAR) * (1 - np.exp(-K1_PAR * LAI))  # Equation 8.27
    R_PAR_FlrCan_up = R_PARGh_ * (1 - np.exp(-K1_PAR * LAI)) * r_FlrPAR * (1 - r_CanPAR) * (1 - np.exp(-K2_PAR * LAI))  # Equation 8.29

    return R_PAR_SunCan_dn + R_PAR_FlrCan_up


def R_NIR_SunCan(h_GlobAir, a_CanNIR, h_GlobNIR, i_Glob):
    """The NIR absorbed by the canopy
    Equation 8.33

    :param h_GlobAir: the ratio of the global radiation
    :param a_CanNIR: NIR absorption coefficient of the canopy
    :param h_GlobNIR: the ratio between NIR and the global radiation
    :param i_Glob: the outside global radiation
    :return: The NIR absorbed by the canopy [W*(m^-2)]
    """
    return (1 - h_GlobAir) * a_CanNIR * h_GlobNIR * i_Glob


# 8.6.3 Convection and conduction
def H(HEC, T1, T2) -> float:
    """Convective and conductive heat fluxes
    Equation 8.40

    :param float HEC: the heat exchange coefficient between object 1 and 2
    :param float T1: the temperature of object 1
    :param float T2: the temperature of object 2
    :return: The heat flow from object 1 to object 2 [W*(m^-2)]
    """
    return HEC * (T1 - T2)

