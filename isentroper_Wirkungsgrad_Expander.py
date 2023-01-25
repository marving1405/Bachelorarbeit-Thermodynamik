def isentroper_Wirkungsgrad(m,n):
    """
        berechnet den isentropen Wirkungsgrad des Expanders
        :param m: mass flow in kg/s
        :param n: rotational speed in 1/s
        :return: isentropic efficiency of expander
        """

    eta_Expander = -10 * m + 1
    return eta_Expander