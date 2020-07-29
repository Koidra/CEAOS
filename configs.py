from math import inf
"""
NOTE: In the GreenLight model, there are no whitewash and shadow screen
"""

class Constants(object):
    def __init__(self):
        pass

    class Global:
        alpha_LeafAir = 5
        delta_H = 2.45E6
        sigma = 5.670E-8
        epsilon_Can = 1
        epsilon_Sky = 1
        eta_GlobNIR = 0.5
        eta_GlobPAR = 0.5
        eta_HeatCO2 = 0.057
        eta_HeatVap = 4.43E-8
        eta_mg_ppm = 0.554
        eta_Roof_Thr = 0.9
        rho_Air0 = 1.20
        rho_CanPAR = 0.07
        rho_CanNIR = 0.35
        rho_Steel = 7850
        rho_Water = 1E3
        gamma = 65.8
        omega = 1.99E-7
        cap_Leaf = 1.2E3
        c_evap1 = 4.30
        c_evap2 = 0.54
        c_day_evap3 = 6.1E-7
        c_night_evap3 = 1.1E-11
        c_day_evap4 = 4.3E-6
        c_night_evap4 = 5.2E-6
        c_pAir = 1E3
        c_pSteel = 0.64E3
        c_pWater = 4.18E3
        g = 9.81
        h_So_j = [0.04, 0.08, 0.16, 0.32, 0.64]
        K_1PAR = 0.7
        K_2PAR = 0.7
        K_NIR = 0.27
        K_FIR = 0.94
        M_Air = 28.96
        M_Water = 18
        R = 8.314E3
        R_Can_SP = 5
        r_b = 275
        r_s_min = 82.0
        s_r_s = -1
        s_MV12 = -0.1

    class Greenhouse:
        class Construction:
            eta_GlobAir = 0.1
            psi = 22
            A_Cov = 9.0E4
            A_Flr = 7.8E4
            c_HECin = 1.86
            c_HECout_1 = 2.8
            c_HECout_2 = 1.2
            c_HECout_3 = 1
            h_Air = 4.7
            h_Elevation = 1470
            h_Gh = 5.1
            C_Gh_d = 0.65
            c_leakage = 1E-4
            C_Gh_w = 0.09
            # h_SideRoof
            h_Vent = 0.97

        class Ventilation:
            # eta_ShScrC_d
            # eta_ShScrC_w
            A_Roof = 7.8E3 # 0.1*A_Flr
            A_Side = 0
            sigma_InsScr = 1
            A_Roof_A_Flr = 0.18
            A_Side_A_Flr = 0

        class Roof:
            epsilon_RfFIR = 0.85
            rho_Rf = 2.6E3
            rho_RfNIR = 0.13
            rho_RfPAR = 0.13
            rho_RfFIR = 0.15
            tau_RfNIR = 0.85
            tau_RfPAR = 0.85
            tau_RfFIR = 0
            lambda_Rf = 1.05
            c_p_Rf = 0.84E3
            h_Rf = 4E-3

        class Whitewash:
            epsilon_ShScrPerFIR = 0.9
            rho_ShScrPer = 1E3
            rho_ShScrPerNIR = 0.3
            rho_ShScrPerPAR = 0.3
            rho_ShScrPerFIR = 0
            tau_ShScrPerNIR = 0.6
            tau_ShScrPerPAR = 0.6
            tau_ShScrPerFIR = 0.1
            lambda_ShScrPer = inf
            c_p_ShScrPer = 4.18E3
            h_ShScrPer = 0.2E-3

        class Shadowscreen:
            # No shadowscreen
            rho_ShScrNIR = 0
            rho_ShScrPAR = 0
            rho_ShScrFIR = 0
            tau_ShScrNIR = 1
            tau_ShScrPAR = 1
            tau_ShScrFIR = 1

        class Thermalscreen:
            epsilon_ThScrFIR = 0.44
            rho_ThScr = 0.2E3
            rho_ThScrNIR = 0.7
            rho_ThScrPAR = 0.7
            rho_ThScrFIR = 0.45
            tau_ThScrNIR = 0.25
            tau_ThScrPAR = 0.25
            tau_ThScrFIR = 0.11
            c_pThScr = 1.8E3
            h_ThScr = 0.35E-3
            K_ThScr = 0.25E-3

        class Blackoutscreen:
            tau_BlScrPAR = 0.01
            rho_BlScrPAR = 0.35

        class Floor:
            epsilon_FlrFIR = 1
            rho_Flr = 2300
            rho_FlrNIR = 0.5
            rho_FlrPAR = 0.65
            lambda_Flr = 1.7
            c_pFlr = 0.88E3
            h_Flr = 0.02

        class Soil:
            rho_c_p_So = 1.73E6
            lambda_soil = 0.85
            h_So = [0.04, 0.08, 0.16, 0.32, 0.64]

        class Heating:
            epsilon_Pipe = 0.88
            phi_Pipe_e = 51E-3
            phi_Pipe_i = 47E-3
            l_Pipe = 1.25

        class ActiveClimateControl:
            # eta_Pad
            # phi_Fog
            # phi_Pad
            # phi_VentForced
            phi_ExtCO2 = 4.3E5
            # COP_MechCool
            # HEC_PasAir
            # P_Blow
            # P_Boil
            # P_Geo
            # P_Ind
            # P_MechCool

        class Lamp:
            # No lamps
            thetaLampMax = 0 # Maximum intensity of lamps
            eta_LampPAR = 0  # fraction of lamp input converted to PAR
            eta_LampNIR = 0  # fraction of lamp input converted to NIR

        class Interlight:
            # No lamps
            thetaIntLampMax = 0 # Maximum intensity of lamps
            eta_IntLampPAR = 0  # fraction of lamp input converted to PAR
            eta_IntLampNIR = 0  # fraction of lamp input converted to NIR