
def cap_Can(cap_Leaf, LAI) -> float:
    """The heat capacity of the canopy
    Equation 8.20
    :param float cap_Leaf: the heat capacity of a square meter of canopy
    :param float LAI: the leaf area index
    :return float: [J*K^(-1)*m^(-2)]
    """
    return cap_Leaf * LAI
