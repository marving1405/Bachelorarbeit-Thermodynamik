# simplified alpha calculation from Nusselt correlation for BA Gertenbach

def alpha_inside_tube(re: object, pr: object, lambda_fluid: object, d_i: object) -> object:
    """
    calculates the heat transfer coefficient for inner tube, only turbulent case
    :rtype: object
    :param re: Reynoldsnumber, dimensionless
    :param pr: Prandtlnumber, dimensionless
    :param lambda_fluid: thermal conductivity of fluid, W/(Km)
    :param d_i: inner diameter, m
    :return: heat transfer coefficient, W/(m**2*K)
    """
    if re < 10000:
        print("ERROR: only applicable for turbulent flow")
    nu = 0.023 * re ** (4/5) * pr ** 0.3                # 0.3 for cooling inner fluidtl-
    alpha = nu * lambda_fluid / d_i
    return alpha


def alpha_outside_tube(d_ai, d_aa, lambda_fluid):
    """
    calculates the heat transfer coefficient for outer tube, only laminar flow
    :param d_ai: inner diameter of outer tube, m
    :param d_aa: outer diameter of outer tube, m
    :param lambda_fluid: thermal conductivity of fluid, W/(mK)
    :return: heat transfer coefficient, W/(m**2*K)
    """
    nu = 3.66 + 1.2 * (d_ai / d_aa) ** (-0.8)
    alpha = nu * lambda_fluid / (d_aa - d_ai)
    return alpha
