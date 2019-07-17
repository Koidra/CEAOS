import numpy as np
from abc import ABC, abstractmethod
from configs import Config as C
from data_models import Setpoint
from equations import capacities, heatfluxes, lumpedcoverlayers


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
        tCovPAR = 1  # Temporary value for tCovPAR
        iGlob = 1  # Temporary value for iGlob

        # Calculating tCovPAR, Section 8.4
        _ShScr_ShScrPer_params = {
            'u1': setpoint.U_ShadingScreen,
            'u2': setpoint.U_ShadingScreenPer,
            't1': C.Greenhouse.LumpedCoverLayers.PAR.t_ShScr,
            't2': C.Greenhouse.LumpedCoverLayers.PAR.t_ShScrPer,
            'r1': C.Greenhouse.LumpedCoverLayers.PAR.r_ShScr,
            'r2': C.Greenhouse.LumpedCoverLayers.PAR.r_ShScrPer
        }

        _Rf_ThScr_params = {
            'u1': setpoint.U_Roof,
            'u2': setpoint.U_ThermalScreen,
            't1': C.Greenhouse.LumpedCoverLayers.PAR.t_Rf,
            't2': C.Greenhouse.LumpedCoverLayers.PAR.t_ThScr,
            'r1': C.Greenhouse.LumpedCoverLayers.PAR.r_Rf,
            'r2': C.Greenhouse.LumpedCoverLayers.PAR.r_ThScr
        }

        t_ShScr_ShScrPer = lumpedcoverlayers.transmission_coefficient_with_control(**_ShScr_ShScrPer_params)
        r_ShScr_ShScrPer = lumpedcoverlayers.reflection_coefficient_with_control(**_ShScr_ShScrPer_params)

        t_Rf_ThScr = lumpedcoverlayers.transmission_coefficient_with_control(**_Rf_ThScr_params)
        r_Rf_ThScr = lumpedcoverlayers.reflection_coefficient_with_control(**_Rf_ThScr_params)

        # From 2 combined layers above
        tCovPAR = lumpedcoverlayers.transmission_coefficient(t_ShScr_ShScrPer, t_Rf_ThScr, r_ShScr_ShScrPer, r_Rf_ThScr)

        # Canopy heat capacity
        capCan = capacities.canopy_heat_capacity(C.Global.capLeaf, LAI)

        # The PAR above the canopy
        R_PAR_Gh = heatfluxes.R_PARGh(C.Global.h_GlobAir, tCovPAR, C.Global.h_GlobPAR, iGlob)

        # -> R_PAR_SunCan
        R_PAR_SunCan = heatfluxes.R_PAR_SunCan(R_PARGh_=R_PAR_Gh,
                                               r_CanPAR=C.Global.r_CanPAR,
                                               r_FlrPAR=C.Greenhouse.r_FlrPAR,
                                               K1_PAR=C.Global.K_1PAR,
                                               K2_PAR=C.Global.K_2PAR,
                                               LAI=LAI)
        # -> R_NIR_SunCan
        # TODO: find a_CanNIR (Equation 8.33)


class VerticalFarmClimateModel(IndoorClimateModel):

    def __init__(self, vertical_farm_struct):
        super(VerticalFarmClimateModel, self).__init__(vertical_farm_struct)

    def step(self, crop_obs: np.ndarray, setpoint: np.ndarray):
        raise NotImplementedError
