from typing import NamedTuple


class Setpoints(NamedTuple):
    U_Blow: float
    U_Boil: float
    U_HeatInd: float
    U_HeatGeo: float
    U_MechCool: float
    U_Fog: float
    U_Roof: float
    U_Side: float
    U_VentForced: float
    U_ExtCO2: float
    U_ShScr: float
    U_ThScr: float
    U_Ind: float
    U_Geo: float
    U_BlScr: float
    U_Lamp: float
    U_IntLamp: float


class States(NamedTuple):
    pipe_t: float
    can_t: float
    air_t: float
    internal_cov_t: float
    external_cov_t: float
    thermal_screen_t: float
    above_thermal_screen_t: float
    floor_t: float
    soil_j_t: [float, float, float, float, float]
    CO2_Air: float
    CO2_Top: float
    # Recheck these vars
    LAI: float
    mechcool_t: float
    MC_AirCan: float


class Weather(NamedTuple):
    I_Glob: float
    outdoor_t: float
    sky_t: float
    soil_out_t: float
    CO2_Out: float
    VP_Out: float
    v_Wind: float
