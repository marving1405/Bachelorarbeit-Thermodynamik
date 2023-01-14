# simplified alpha calculation from Nusselt correlation for BA Gertenbach
import ht
import CoolProp.CoolProp as CP
import numpy as np

def alpha_1P_i(p, T, fluid, m_dot, d_i):
    """
    calculates heat transfer coefficient inside tube for 1 phase
    :param p: pressure Pa
    :param T: temperature K
    :param fluid: fluid as string
    :param m_dot: mass flow in kg/s
    :param d_i: inner diameter
    :return: heat transfer coefficient W/(m^2*K)
    """
    vis = CP.PropsSI("VISCOSITY", "P", p, "T", T, fluid)
    rho = CP.PropsSI("D", "P", p, "T", T, fluid)
    Pr = CP.PropsSI("PRANDTL", "P", p, "T", T, fluid)
    lam = CP.PropsSI("L", "P", p, "T", T, fluid)
    Re = 4 * m_dot / (vis * np.pi * d_i)
    Nu = ht.conv_internal.Nu_conv_internal(Re, Pr, Di=d_i)
    alpha = lam * Nu / d_i
    return alpha


def alpha_1P_annulus(p, T, fluid, m_dot, d_ai, d_a):
    """
    calculates heat transfer coefficient inside tube annulus for 1 phase
    :param p: pressure
    :param T: temperature
    :param fluid: fluid as string
    :param m_dot: massflow in kg/s
    :param d_ai: inner diameter of outer tube
    :param d_a: outer diameter of outer tube
    :return: heat transfer coefficient in W/(m^2*K)
    """
    vis = CP.PropsSI("VISCOSITY", "P", p, "T", T, fluid)
    rho = CP.PropsSI("D", "P", p, "T", T, fluid)
    Pr = CP.PropsSI("PRANDTL", "P", p, "T", T, fluid)
    lam = CP.PropsSI("L", "P", p, "T", T, fluid)
    Re = 4 * m_dot * (d_a - d_ai) / (vis * np.pi * (d_a ** 2 - d_ai ** 2))
    if Re > 2300:
        Nu = 0.023 * Re ** 0.8 * Pr ** 0.4
    else:
        Nu = 3.66 + 1.2 * (d_ai / d_a) ** -0.8
    alpha = lam * Nu / (d_a - d_ai)
    return alpha

def alpha_boiling(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te):
    """
    calculates alpha in wet steam region, documentation available on https://ht.readthedocs.io/en/latest/ht.boiling_flow.html
    :param m: Massenstrom
    :param x: Dampfgehalt im Rohrintervall
    :param D: Rohrdurchmesser
    :param rhol: Dichte Flüssigkeit
    :param rhog: Dichte Gas
    :param mul: Viskosität Flüssigkeit
    :param mug: Viskosität Gas
    :param kl: Wärmeleitfähigkeit lambda
    :param Cpl: Wärmekapazität Flüssigkeit
    :param Hvap: Verdampfungsenthalpie Flüssigkeit
    :param sigma: Oberflächenspannung Flüssigkeit
    :param dPsat: Differenz des Sättigungsdruckes der Flüssigkeit am Einlass und Auslass
    :param Te: Übertemperatur der Wand
    :return:
    """
    alpha = ht.boiling_flow.Chen_Bennett(m, x, D, rhol, rhog, mul, mug, kl, Cpl, Hvap, sigma, dPsat, Te)
    return alpha


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
