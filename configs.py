
class GHCommonConstants:
    h_GlobAir = 0.1
    """
    η_GlobAir, table 8.1
    The ratio of the global radiation which is absorbed by the greenhouse construction elements
    """

    h_GlobPAR = 0.5
    """
    η_GlobPAR, table 8.1
    Ratio between PAR and the outside global radiation
    """

    h_GlobNIR = 0.5
    """
    ηGlob_ NIR
    Ratio between NIR and the outside global radiation 
    """

    r_CanPAR = 0.07
    """
    ρ_CanPAR, table 8.1
    The PAR reflection coefficient
    """

    r_CanNIR = 0.35
    """
    ρCanNIR, table 8.1
    The NIR reflection coefficient of the top of the canopy
    """

    K_1PAR = 0.7
    """
    K_1PAR, table 8.1
    PAR extinction coefficient of the canopy
    """

    K_2PAR = 0.7
    """
    K_2PAR, table 8.1
    PAR extinction coefficient of the canopy when PAR is reflected from the floor 
    """

    K_NIR = 0.27
    """
    K_NIR, table 8.1
    Extinction coefficient of the canopy for NIR
    """

    a_LeafAir = 5
    """
    αLeafAir [W.m^(-2).K^(-1)], table 8.1
    Convective heat exchange coefficient from the canopy leaf to the greenhouse air 
    """

    d_H = 2.45 * 10E6
    """
    ∆H, table 8.1
    Latent heat of evaporation
    """

    c_pAir = 10E3
    """
    c_pAir, table 8.1
    Specific heat capacity of the air
    """

    r_FlrPAR = 0.65  # Texas
    """
    ρ_FlrPAR, table 8.1
    PAR reflection coefficient of the floor
    """

    r_FlrNIR = 0.5  # Texas
    """
    ρ_FlrNIR, table 8.1
    NIR reflection coefficient of the floor
    """


class LumpedCoverConstants:
    """
    This class is based on section 8.4 and table 8.2 (in section 8.13).
    The model contains four cover layers, i.e.
    a movable outdoor shading screen (ShScr),
    a semi-permanent shading screen (ShScrPer),
    the greenhouse roof (Rf) and
    a movable indoor thermal screen (ThScr).

    Default value is of a Texas greenhouse.
    Shading Screen's constants are not mentioned in the Thesis, so we set reflection coefficient to 0
    and the transmission coefficient to 1.

    Due to the complexity of parameters,
    this class hold temporarily value of Lumped cover layers.
    All of these configs will be exported to YAML file.


    r_: The reflection coefficient of the layer
    t_: The transmission coefficient of the layer
    """
    class PAR:
        # Shading Screen
        r_ShScr = 0
        t_ShScr = 1

        # Semi-permanent shading screen
        r_ShScrPer = 0.3
        t_ShScrPer = 0.6

        # Roof layer
        r_Rf = 0.13
        t_Rf = 0.85

        # Thermal Screen
        r_ThScr = 0.7
        t_ThScr = 0.25

    class NIR:
        # Shading Screen
        r_ShScr = 0
        t_ShScr = 1

        # Semi-permanent shading screen
        r_ShScrPer = 0.3
        t_ShScrPer = 0.6

        # Roof layer
        r_Rf = 0.13
        t_Rf = 0.85

        # Thermal Screen
        r_ThScr = 0.7
        t_ThScr = 0.25

    class FIR:
        # Shading Screen
        r_ShScr = 0
        t_ShScr = 1

        # Semi-permanent shading screen
        r_ShScrPer = 0
        t_ShScrPer = 0.1

        # Roof layer
        r_Rf = 0.15
        t_Rf = 0

        # Thermal Screen
        r_ThScr = 0.45
        t_ThScr = 0.11
