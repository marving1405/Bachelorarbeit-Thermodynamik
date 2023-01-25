import numpy as np

def lambda_Oel(T):
    """
        calculates Thermal Conductivity of shell heat transfer oil s2
        :param T: temperature of oil at input of high temperature tank
        """
    T_umrechnung = T - 273.15
    lambda_oel = -7.22E-5 * T_umrechnung + 0.136
    return lambda_oel


def cp_Oel(T):
    """
        calculates Thermal Conductivity of shell heat transfer oil s2
        :param T: temperature of oil at input of high temperature tank
        """
    T_umrechnung = T - 273.15
    cp_oel = 0.004 * T_umrechnung + 1.809
    return cp_oel