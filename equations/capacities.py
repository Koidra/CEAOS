
def canopy_heat_capacity(capLeaf: float, LAI: float) -> float:
    """The heat capacity of the canopy
    Equation 8.20

    :param capLeaf: the heat capacity of a square meter of canopy
    :param LAI: the leaf area index
    :return: [J K^(-1) m^(-2)]
    """
    return capLeaf * LAI
