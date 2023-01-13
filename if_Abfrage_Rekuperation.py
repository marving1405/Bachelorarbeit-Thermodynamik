import json, CoolProp.CoolProp as CP
from calculate_alpha_aw import alpha_1P_i, alpha_inside_tube, alpha_outside_tube
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
from CoolProp.CoolProp import PropsSI
import numpy as np

fluid = "REFPROP::PROPANE"

def abfrage(T4, p4, m_ORC, A_quer, d_i, d_ai, d_aa, h4, A_i, A_a):
    abfrage = input("Rekuperator verwenden? True or False?")
    if abfrage == "True":
        # Rekuperator (ersetzt WÜ1 und K1)

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

    if abfrage == "False":
        # Kondensator 1, ÜD -> SF

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


