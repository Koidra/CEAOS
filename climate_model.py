import numpy as np
from abc import ABC, abstractmethod
from configs import Config as C
from data_models import Setpoint
from equations import capacities, heatfluxes, lumpedcoverlayers
from utils import params_except_keys


class IndoorClimateModel(ABC):
    """ An indoor climate model
    Takes a climate control setpoint vector and returns a climate observation vector as output
    """

    # Greenhouse or Vertical farm design constants
    # E.g.: Lamps Par, Light Transmission Ratio, ...
    greenhouse_design = None

    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def step(self, crop_obs: np.ndarray, setpoint: Setpoint):
        pass


class GreenhouseClimateModel(IndoorClimateModel):

    def __init__(self, greenhouse_config_file):
        super(GreenhouseClimateModel, self).__init__(greenhouse_config_file)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError

    def __canopy_temperature(self, setpoint: Setpoint):
        LAI = 1  # Temporary value of Leaf Area Index
        iGlob = 1  # Temporary value for iGlob

        # ***************************************************
        # Calculating tCovPAR, Section 8.4
        # ***************************************************

        _ShScr_ShScrPer_PAR_params = {
            'u1': setpoint.U_ShadingScreen,
            'u2': setpoint.U_ShadingScreenPer,
            't1': C.Greenhouse.LumpedCoverLayers.PAR.t_ShScr,
            't2': C.Greenhouse.LumpedCoverLayers.PAR.t_ShScrPer,
            'r1': C.Greenhouse.LumpedCoverLayers.PAR.r_ShScr,
            'r2': C.Greenhouse.LumpedCoverLayers.PAR.r_ShScrPer
        }

        _Rf_ThScr_PAR_params = {
            'u1': setpoint.U_Roof,
            'u2': setpoint.U_ThermalScreen,
            't1': C.Greenhouse.LumpedCoverLayers.PAR.t_Rf,
            't2': C.Greenhouse.LumpedCoverLayers.PAR.t_ThScr,
            'r1': C.Greenhouse.LumpedCoverLayers.PAR.r_Rf,
            'r2': C.Greenhouse.LumpedCoverLayers.PAR.r_ThScr
        }

        t_PAR_ShScr_ShScrPer = lumpedcoverlayers.transmission_coefficient(**_ShScr_ShScrPer_PAR_params)
        r_PAR_ShScr_ShScrPer = lumpedcoverlayers.reflection_coefficient(**params_except_keys(_ShScr_ShScrPer_PAR_params, except_keys='t2'))

        t_PAR_Rf_ThScr = lumpedcoverlayers.transmission_coefficient(**_Rf_ThScr_PAR_params)
        r_PAR_Rf_ThScr = lumpedcoverlayers.reflection_coefficient(**params_except_keys(_Rf_ThScr_PAR_params, except_keys='t2'))

        # From 2 combined layers above
        tCovPAR = lumpedcoverlayers.transmission_coefficient(t_PAR_ShScr_ShScrPer, t_PAR_Rf_ThScr, r_PAR_ShScr_ShScrPer, r_PAR_Rf_ThScr)

        # ***************************************************
        # Calculate r_Cov_NIR, Using in equation 8.30
        # ***************************************************

        _ShScr_ShScrPer_NIR_params = {
            'u1': setpoint.U_ShadingScreen,
            'u2': setpoint.U_ShadingScreenPer,
            't1': C.Greenhouse.LumpedCoverLayers.NIR.t_ShScr,
            't2': C.Greenhouse.LumpedCoverLayers.NIR.t_ShScrPer,
            'r1': C.Greenhouse.LumpedCoverLayers.NIR.r_ShScr,
            'r2': C.Greenhouse.LumpedCoverLayers.NIR.r_ShScrPer
        }

        _Rf_ThScr_NIR_params = {
            'u1': setpoint.U_Roof,
            'u2': setpoint.U_ThermalScreen,
            't1': C.Greenhouse.LumpedCoverLayers.NIR.t_Rf,
            't2': C.Greenhouse.LumpedCoverLayers.NIR.t_ThScr,
            'r1': C.Greenhouse.LumpedCoverLayers.NIR.r_Rf,
            'r2': C.Greenhouse.LumpedCoverLayers.NIR.r_ThScr
        }

        t_NIR_ShScr_ShScrPer = lumpedcoverlayers.transmission_coefficient(**_ShScr_ShScrPer_NIR_params)
        r_NIR_ShScr_ShScrPer = lumpedcoverlayers.reflection_coefficient(**params_except_keys(_ShScr_ShScrPer_NIR_params, except_keys='t2'))
        r_NIR_Rf_ThScr = lumpedcoverlayers.reflection_coefficient(**params_except_keys(_Rf_ThScr_NIR_params, except_keys='t2'))

        # From 2 combined layers above
        rCovNIR = lumpedcoverlayers.reflection_coefficient(t_NIR_ShScr_ShScrPer, r_NIR_ShScr_ShScrPer, r_NIR_Rf_ThScr)
        tCovNIR_virtual = 1 - rCovNIR  # Equation 8.30

        # ***************************************************
        # Calculate NIR absorption coefficient of the canopy, Using in equation 8.33
        # ***************************************************

        tCanNIR_virtual = np.exp(-C.Global.K_NIR*LAI)  # Equation 8.31
        rCanNIR_virtual = C.Global.r_CanNIR*(1-tCanNIR_virtual)  # Equation 8.32

        tMultiLayerCanNIR = lumpedcoverlayers.transmission_coefficient(tCovNIR_virtual, tCanNIR_virtual, rCovNIR, rCanNIR_virtual, u2=setpoint.U_Canopy)
        rMultiLayerCanNIR = lumpedcoverlayers.reflection_coefficient(tCovNIR_virtual, rCovNIR, rCanNIR_virtual, u2=setpoint.U_Canopy)

        aCanNIR = lumpedcoverlayers.absorption_coefficient(tMultiLayerCanNIR, rMultiLayerCanNIR)  # in equation 8.33

        # ---> R_NIR_SunCan
        R_NIR_SunCan = heatfluxes.R_NIR_SunCan(h_GlobAir=C.Global.h_GlobAir,
                                               a_CanNIR=aCanNIR,
                                               h_GlobNIR=C.Global.h_GlobNIR,
                                               i_Glob=iGlob)

        # ***************************************************
        # Calculate NIR absorption coefficient of the floor, Using in equation 8.34
        # ***************************************************
        tFlrNIR_virtual = 1 - C.Greenhouse.r_FlrNIR  # Equation 8.30

        tMultiLayerFlrNIR = lumpedcoverlayers.transmission_coefficient(tCovNIR_virtual, tFlrNIR_virtual, rCovNIR, C.Greenhouse.r_FlrNIR, u2=setpoint.U_Floor)
        rMultiLayerFlrNIR = lumpedcoverlayers.reflection_coefficient(tCovNIR_virtual, rCovNIR, C.Greenhouse.r_FlrNIR, u2=setpoint.U_Floor)

        aFlrNIR = lumpedcoverlayers.absorption_coefficient(tMultiLayerFlrNIR, rMultiLayerFlrNIR)  # Using in equation 8.34

        # ***************************************************
        # Canopy heat capacity
        # ***************************************************
        capCan = capacities.canopy_heat_capacity(C.Global.capLeaf, LAI)

        # ***************************************************
        # The PAR above the canopy
        # ***************************************************
        R_PAR_Gh = heatfluxes.R_PARGh(C.Global.h_GlobAir, tCovPAR, C.Global.h_GlobPAR, iGlob)

        # -> R_PAR_SunCan
        R_PAR_SunCan = heatfluxes.R_PAR_SunCan(R_PARGh_=R_PAR_Gh,
                                               r_CanPAR=C.Global.r_CanPAR,
                                               r_FlrPAR=C.Greenhouse.r_FlrPAR,
                                               K1_PAR=C.Global.K_1PAR,
                                               K2_PAR=C.Global.K_2PAR,
                                               LAI=LAI)
        # -> R_NIR_SunCan
        # TODO: validate a_CanNIR (Equation 8.33)


class VerticalFarmClimateModel(IndoorClimateModel):

    def __init__(self, vertical_farm_struct):
        super(VerticalFarmClimateModel, self).__init__(vertical_farm_struct)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError
