import numpy as np
from .constanst import LumpedCoverConstants as LCC
from .constanst import GHCommonConstants as C


def __multiply(**kwargs):
    return np.prod(np.array(list(kwargs.values())))


class LumpedCoverLayers:
    """
    8.4
    The model contains four cover layers, i.e.
    a movable outdoor shading screen (ShScr),
    a semi-permanent shading screen (ShScrPer),
    the greenhouse roof (Rf) and
    a movable indoor thermal screen (ThScr).
    """

    __consts = None
    __U_ShScr = 1
    __U_ShScrPer = 1
    __U_Rf = 1
    __U_ThScr = 1

    def __init__(self, U_ShScr, U_ShScrPer, U_Rf, U_ThScr):
        """
        :param float U_ShScr: The control of the movable shading screen
        :param float U_ShScrPer: The control of the semi-permanent shading screen
        :param float U_Rf: The control of the roof
        :param float U_ThScr: The control of the thermal screen
        """

        self.__U_ShScr = U_ShScr
        self.__U_ShScrPer = U_ShScrPer
        self.__U_Rf = U_Rf
        self.__U_ThScr = U_ThScr

    def __parse_consts(self, target_radiation):
        if target_radiation == 'PAR':
            return LCC.PAR
        if target_radiation == 'NIR':
            return LCC.NIR
        if target_radiation == 'FIR':
            return LCC.FIR

    def _t_ShScr_ShScrPer(self, target='PAR'):
        """
        Equation 8.16
        """
        consts = self.__parse_consts(target)
        return ((1 - self.__U_ShScr * (1 - consts.t_ShScr)) * (
                1 - self.__U_ShScrPer * (1 - consts.t_ShScrPer))) / \
               (1 - self.__U_ShScr * consts.r_ShScr * self.__U_ShScrPer * consts.r_ShScrPer)

    def _r_ShScr_ShScrPer(self, target='PAR'):
        """
        Equation 8.17
        """
        consts = self.__parse_consts(target)
        return self.__U_ShScr * consts.r_ShScr + \
               (((1 - self.__U_ShScr * (
                       1 - consts.t_ShScr)) ** 2) * self.__U_ShScrPer * consts.r_ShScrPer) / \
               (1 - self.__U_ShScr * consts.r_ShScr * self.__U_ShScrPer * consts.r_ShScrPer)

    def _t_Rf_ThScr(self, target='PAR'):
        """
        Page 208: Secondly,...
        """
        consts = self.__parse_consts(target)
        return ((1 - self.__U_Rf * (1 - consts.t_Rf)) * (1 - self.__U_ThScr * (1 - consts.t_ThScr))) / \
               (1 - self.__U_Rf * consts.r_Rf * self.__U_ThScr * consts.r_ThScr)

    def _r_Rf_ThScr(self, target='PAR'):
        """
        Page 208: Secondly,...
        """
        consts = self.__parse_consts(target)
        return self.__U_Rf * consts.r_Rf + \
               (((1 - self.__U_Rf * (1 - consts.t_Rf)) ** 2) * self.__U_ThScr * consts.r_ThScr) / \
               (1 - self.__U_Rf * consts.r_Rf * self.__U_ThScr * consts.r_ThScr)

    def t_Cover(self, target='PAR'):
        """
        The transmission coefficient
        Equation 8.14
        """

        __t_ShScr_ShScrPer = self._t_ShScr_ShScrPer(target)
        __r_ShScr_ShScrPer = self._r_ShScr_ShScrPer(target)
        __t_Rf_ThScr = self._t_ShScr_ShScrPer(target)
        __r_Rf_ThScr = self._r_Rf_ThScr(target)

        return (__t_ShScr_ShScrPer * __t_Rf_ThScr) / (1 - __r_ShScr_ShScrPer * __r_Rf_ThScr)

    def r_Cover(self, target='PAR'):
        """
        the reflection coefficient
        Equation 8.15
        """

        __t_ShScr_ShScrPer = self._t_ShScr_ShScrPer(target)
        __r_ShScr_ShScrPer = self._r_ShScr_ShScrPer(target)
        __t_Rf_ThScr = self._t_ShScr_ShScrPer(target)
        __r_Rf_ThScr = self._r_Rf_ThScr(target)

        return __r_ShScr_ShScrPer + ((__t_ShScr_ShScrPer ** 2) * __r_Rf_ThScr) / (1 - __r_ShScr_ShScrPer * __r_Rf_ThScr)


class HeatFluxes:
    """
    8.6
    """
    setpoint = {}

    def __init__(self, setpoint):
        self.setpoint = setpoint
        self._lumped_cover_layers = LumpedCoverLayers(U_ShScr=self.setpoint['U_ShScr'],
                                                      U_ShScrPer=self.setpoint['U_ShScrPer'],
                                                      U_Rf=self.setpoint['U_Rf'],
                                                      U_ThScr=self.setpoint['U_ThScr'])

    # 8.3.1 Global, PAR and NIR heat fluxes

    def R_PAR_Gh(self, i_Glob):
        """
        The PAR above the canopy
        Equation 8.28
        :param i_Glob:
        :return: R_PAR_Gh [W.m^-2] - equation 8.28
        """
        t_CovPAR = self._lumped_cover_layers.t_Cover('PAR')
        return (1 - C.h_GlobAir) * t_CovPAR * C.h_GlobPAR * i_Glob

    def R_NIR_Gh(self, i_Glob):
        """
        The NIR above the canopy
        Based on equation 8.28
        :param i_Glob:
        """
        t_CovNIR = self._lumped_cover_layers.t_Cover('NIR')
        return (1 - C.h_GlobAir) * t_CovNIR * C.h_GlobNIR * i_Glob

    def R_PAR_SunCan(self, i_Glob, LAI):
        _R_PAR_Gh = self.R_PAR_Gh(i_Glob)
        _R_PAR_SunCan_down = _R_PAR_Gh * (1 - C.r_CanPAR) * (1 - np.exp(-C.K_1PAR * LAI))
        _R_PAR_FlrCan_up = _R_PAR_Gh * np.exp(-C.K_1PAR * LAI) * C.r_FlrPAR * (1 - C.r_CanPAR) * (
                1 - np.exp(-C.K_2PAR * LAI))

        return _R_PAR_SunCan_down + _R_PAR_FlrCan_up

    def R_NIR_SunCan(self, i_Glob, LAI):
        _R_NIR_Gh = self.R_NIR_Gh(i_Glob)
        _R_NIR_SunCan_down = _R_NIR_Gh * (1 - C.r_CanNIR) * (1 - np.exp(-C.K_NIR * LAI))
        _R_NIR_FlrCan_up = _R_NIR_Gh * np.exp(-C.K_NIR * LAI) * C.r_FlrNIR * (1 - C.r_CanNIR) * (
                1 - np.exp(-C.K_NIR * LAI))

        return _R_NIR_SunCan_down + _R_NIR_FlrCan_up

    # 8.6.3 Convective and conductive heat fluxes
    def HEC_CanAir(self, LAI):
        return 2 * C.a_LeafAir * LAI

    def H_CanAir(self, LAI, T_Can, T_Air):
        """
        Equation 8.40
        The sensible heat exchange between canopy and greenhouse air
        :param LAI: Leaf Area Index
        :param T_Can: The temperature of Canopy
        :param T_Air: The temperature of Air
        """
        return self.HEC_CanAir(LAI)*(T_Can - T_Air)

    def delta_TCan(self, LAI, T_Can, T_Air):
        pass

    def delta_TAir(self, LAI, T_Can, T_Air):
        pass

    # 8.6.4 Latent heat fluxes
    def LatenHeatFluxes(self):
        pass

    # 8.9 Canopy transpiration
    def VEC_CanAir(self, LAI, r_Air):
        pass


class Calculation:
    setpoint = {}

    def __init__(self, **kwargs):
        self.setpoint = kwargs
        self._lumped_cover_layers = LumpedCoverLayers(U_ShScr=self.setpoint['U_ShScr'],
                                                      U_ShScrPer=self.setpoint['U_ShScrPer'],
                                                      U_Rf=self.setpoint['U_Rf'],
                                                      U_ThScr=self.setpoint['U_ThScr'])

    def R_PAR_Gh(self, i_Glob):
        """
        The PAR above the canopy
        Equation 8.28
        :param i_Glob:
        :return: R_PAR_Gh [W.m^-2] - equation 8.28
        """
        t_CovPAR = self._lumped_cover_layers.t_Cover('PAR')
        return (1 - C.h_GlobAir) * t_CovPAR * C.h_GlobPAR * i_Glob

    def R_NIR_Gh(self, i_Glob):
        """
        The NIR above the canopy
        Based on equation 8.28
        :param i_Glob:
        """
        t_CovNIR = self._lumped_cover_layers.t_Cover('NIR')
        return (1 - C.h_GlobAir) * t_CovNIR * C.h_GlobNIR * i_Glob

    def R_PAR_SunCan(self, i_Glob, LAI):
        _R_PAR_Gh = self.R_PAR_Gh(i_Glob)
        _R_PAR_SunCan_down = _R_PAR_Gh * (1 - C.r_CanPAR) * (1 - np.exp(-C.K_1PAR * LAI))
        _R_PAR_FlrCan_up = _R_PAR_Gh * np.exp(-C.K_1PAR * LAI) * C.r_FlrPAR * (1 - C.r_CanPAR) * (
                1 - np.exp(-C.K_2PAR * LAI))

        return _R_PAR_SunCan_down + _R_PAR_FlrCan_up

    def R_NIR_SunCan(self, i_Glob, LAI):
        _R_NIR_Gh = self.R_NIR_Gh(i_Glob)
        _R_NIR_SunCan_down = _R_NIR_Gh * (1 - C.r_CanNIR) * (1 - np.exp(-C.K_NIR * LAI))
        _R_NIR_FlrCan_up = _R_NIR_Gh * np.exp(-C.K_NIR * LAI) * C.r_FlrNIR * (1 - C.r_CanNIR) * (
                1 - np.exp(-C.K_NIR * LAI))

        return _R_NIR_SunCan_down + _R_NIR_FlrCan_up
