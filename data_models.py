from typing import NamedTuple


class Setpoints(NamedTuple):
    U_Blow: float
    U_Boil: float
    U_HeatInd: float
    U_HeatGeo: float
    U_Pad: float
    U_MechCool: float
    U_Fog: float
    U_Roof: float
    U_Side: float
    U_VentForced: float
    U_ExtCO2: float
    U_ShScr: float
    U_ShScrPer: float
    U_ThScr: float
    U_Ind: float
    U_Geo: float
    U_BlScr: float
    U_Lamp: float
    U_IntLamp: float

class Inputs(NamedTuple):
    outdoor_t: float
    sky_t: float
    mechcool_t: float
    pipe_t: float
    can_t: float
    air_t: float
    internal_cov_t: float
    external_cov_t: float
    thermal_screen_t: float
    above_thermal_screen_t: float
    floor_t: float
    soil_j_t: [float, float, float, float, float]
    soil_out_t: float
    CO2_Out: float
    CO2_Air: float
    CO2_Top: float
    VP_Out: float
    I_Glob: float
    v_Wind: float
    LAI: float
    MC_AirCan: float
    # soil_mean_t: float
    # t: float
    # a_0: float
    # beta: float
