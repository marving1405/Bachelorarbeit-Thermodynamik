import json, CoolProp.CoolProp as CP
from calculate_alpha_aw import alpha_1P_i, alpha_inside_tube, alpha_outside_tube
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
from CoolProp.CoolProp import PropsSI
import numpy as np
from Solver_fsolve import solveT2_Rekuperator


fluid = "REFPROP::PROPANE"

def abfrage(T4, p4, m_ORC, A_quer, d_i, d_ai, d_aa, h4, A_i, A_a, p2, T2):
    abfrage = input("Rekuperator verwenden? True or False?")
    if abfrage == "True":
    # Rekuperator (ersetzt WÜ1 und K1)


    if abfrage == "False":
        # Kondensator 1 und WÜ1, ÜD -> SF

        # WÜ1
        T2_siedend = PropsSI('T', 'P', p2, 'Q', 0, fluid)
        cp_fluid_1 = PropsSI('C', 'T', T2, 'P', p2, fluid)
        Q_zu1 = m_ORC * cp_fluid_1 * (T2_siedend - T2)
        A_quer = np.pi * (d_i / 2) ** 2  # m2

        rho_1 = PropsSI('D', 'T', T2, 'P', p2, fluid)  # kg/m3
        v_1 = m_ORC / rho_1
        c_1 = v_1 / A_quer
        viscosity_1 = PropsSI('VISCOSITY', 'T', T2, 'P', p2, fluid)

        # re_1 = rho_1 * c_1 * d_i / viscosity_1
        lambda_fluid_1 = PropsSI('CONDUCTIVITY', 'T', T2, 'P', p2, fluid)
        pr_1 = PropsSI('PRANDTL', 'T', T2, 'P', p2, fluid)
        # alpha_i_1 = alpha_inside_tube(re_1,pr_1,lambda_fluid_1, d_i)
        alpha_a_1 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_1)

        from calculate_alpha_aw import alpha_1P_i

        alpha_i_1 = alpha_1P_i(p2, T2, fluid, m_ORC, d_i)

        Tlow_H = 60 + 273.15  # K
        Tlow_L = 20 + 273.15  # K
        dTA_1 = T2_siedend - T2
        dTB_1 = Tlow_H - Tlow_L
        l1 = Q_zu1 / (np.pi * d_i * alpha_i_1 * ((dTA_1) - (dTB_1) / np.log(dTA_1 / dTB_1)))

        A_i = 2 * np.pi * d_i / 2 * l1 / 3
        A_a = 2 * np.pi * d_ai / 2 * l1 / 3
        R_konv_innen1 = 1 / A_i * alpha_i_1
        R_konv_aussen1 = 1 / A_a * alpha_a_1
        R_waermeleitung1 = np.log(d_aa / d_ai) / (2 * np.pi * l1 * lambda_fluid_1)
        R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1

        # K1
        rho_k1 = PropsSI('D', 'T', T4, 'P', p4, fluid)  # kg/m3
        v_k1 = m_ORC / rho_k1
        c_k1 = v_k1 / A_quer
        viscosity_k1 = PropsSI('VISCOSITY', 'T', T4, 'P', p4, fluid)

        re_k1 = rho_k1 * c_k1 * d_i / viscosity_k1
        lambda_fluid_k1 = PropsSI('CONDUCTIVITY', 'T', T4, 'P', p4, fluid)
        pr_k1 = PropsSI('PRANDTL', 'T', T4, 'P', p4, fluid)

        alpha_i_k1 = alpha_inside_tube(re_k1, pr_k1, lambda_fluid_k1, d_i)
        alpha_a_k1 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k1)

        h4_siedend = PropsSI('H', 'P', p4, 'Q', 0, fluid)
        T4_siedend = PropsSI('T', 'P', p4, 'H', h4_siedend, fluid)
        Q_ab1 = m_ORC * (h4 - h4_siedend)

        dTA_k1 = T4 - T4_siedend
        dTB_k1 = 40  # Kühlwasser etc.
        l_k1 = Q_ab1 / (np.pi * d_i * alpha_i_k1 * ((dTA_k1) - (dTB_k1) / np.log(dTA_k1 / dTB_k1)))

        A_i_k1 = 2 * np.pi * d_i / 2 * l_k1 / 2
        A_a_k1 = 2 * np.pi * d_ai / 2 * l_k1
        R_konv_innen = 1 / A_i * alpha_i_k1
        R_konv_aussen = 1 / A_a * alpha_a_k1
        R_waermeleitung = np.log(d_aa / d_ai) / (2 * np.pi * l_k1 * lambda_fluid_k1)
        alpha_R = alpha_1P_i(p4, T4, fluid, m_ORC, d_i)


