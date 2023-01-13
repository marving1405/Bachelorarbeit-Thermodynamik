"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""
import json, CoolProp.CoolProp as CP

CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from calculate_alpha_aw import alpha_inside_tube, alpha_outside_tube

fluid = "REFPROP::PROPANE"
m_ORC = 10E-3

"""
Pumpe: Zustand 1 so gewählt, dass Fluid bei 1 bar und unterkühlt vorliegt

"""
p1 = 100000  # Pa
T1 = 229  # Kelvin
etaP = 0.9  # wird hier als konstant angesehen

# Berechnung der zu verrichtenden Pumpenarbeit
# Eintrittszustand bei 1 bar unterkühlte Flüssigkeit

h1 = PropsSI('H', 'T', T1, 'P', p1, fluid)
p2 = 2000000  # Pa
v1 = 1 / (PropsSI('D', 'P', p1, 'T', T1, fluid))  # m3/kg

w_p = (v1 * (p2 - p1))  # J
P_p = (w_p * m_ORC) / (etaP)
h2 = w_p + h1
T2 = PropsSI('T', 'H', h2, 'P', p2, fluid)

""" Verdampfer
    Betrachtet wird ein Doppelrohrwärmeübertrager
    Als Rohrdurchmesser werden Innen 10mm und außen 14mm festgelegt
    Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
    """

d_i = 10E-3  # m
d_ai = 14E-3
d_aa = 16E-3
A_quer = np.pi * (d_i / 2) ** 2  # m2

from Test_fsolve import solveT3


# TODO Implementieren Massenstromverhältnis
# for x in np.arange(0.1,10,0.1):
m_oel = 4 * m_ORC

cp_oel = 1.9  # kJ/kg*K
Tmittel_H = 100 + 273.15  # K

T4 = 238.9
h4 = 537000
# Rekuperator
# Kühlen ÜD -> SF
from calculate_alpha_aw import alpha_1P_i
from scipy.optimize import fsolve
from Test_fsolve import solveT2_Rekuperator
p4 = p1
lR = 2  # m
h4_sattdampf = PropsSI('H', 'P', p4, 'Q', 1, fluid)
T4_sattdampf = PropsSI('T', 'P', p4, 'H', h4_sattdampf, fluid)
dTA = float(T4 - T4_sattdampf)
alpha_R = float(alpha_1P_i(p4, T4, fluid, m_ORC, d_i))
Q_abR = float(m_ORC * (h4 - h4_sattdampf))
Q_zu1 = Q_abR
T2_R = fsolve(solveT2_Rekuperator, 230., args=(Q_abR, d_i, lR, alpha_R, T2, dTA))
h2_R = PropsSI('H', 'P', p2, 'T', T2_R, fluid)

'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''
#T2_siedend = PropsSI('T', 'P', p2, 'Q', 0, fluid)
lambda_fluid_2 = PropsSI('CONDUCTIVITY', 'T', T2_R, 'Q', 0, fluid)
alpha_i_zweiphasig = 600
alpha_a_zweiphasig = 400
# TODO andere alphas aus VDI Waermeatlas ok?

h2_sattdampf = PropsSI('H', 'P', p2, 'T', T2, fluid)

Q_zu2 = m_ORC * (h2_sattdampf - h2_R)
Tmittel_L = Tmittel_H - ((Q_zu2 / 1000) / (m_oel * cp_oel))

l2 = Q_zu2 / (np.pi * d_i * alpha_i_zweiphasig * (Tmittel_H - Tmittel_L))
A_i_2 = 2 * np.pi * d_i / 2 * l2 / 3
A_a_2 = 2 * np.pi * d_ai / 2 * l2 / 3

R_konv_innen2 = 1 / A_i_2 * alpha_i_zweiphasig
R_konv_aussen2 = 1 / A_a_2 * alpha_a_zweiphasig
R_waermeleitung2 = np.log(d_aa / d_ai) / (2 * np.pi * l2 * lambda_fluid_2)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2

'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
'''

Thoch_H = 180 + 273.15  # K
Thoch_L = 120 + 273.15  # K
l3 = 5  # m
Q_zu3 = m_oel * cp_oel * (Thoch_H - Thoch_L) * 1000
# h3 = PropsSI('H','T',T3,'P',p2,fluid)
T2_sattdampf = PropsSI('T', 'P', p2, 'Q', 1, fluid)

rho_3 = PropsSI('D', 'T', T2_sattdampf, 'Q', 1, fluid)  # kg/m3
v_3 = m_ORC / rho_3
c_3 = v_3 / A_quer

viscosity_3 = PropsSI('VISCOSITY', 'T', T2_sattdampf, 'Q', 1, fluid)

re_3 = rho_3 * c_3 * d_i / viscosity_3
lambda_fluid_3 = PropsSI('CONDUCTIVITY', 'T', T2_sattdampf, 'Q', 1, fluid)
pr_3 = PropsSI('PRANDTL', 'T', T2_sattdampf, 'Q', 1, fluid)
alpha_i_3 = alpha_inside_tube(re_3, pr_3, lambda_fluid_3, d_i)  # alpha_1P_i(p2, T2_sattdampf, fluid, m_ORC, d_i)
alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_3)

A_i_3 = 2 * np.pi * d_i / 2 * l3 / 3
A_a_3 = 2 * np.pi * d_ai / 2 * l3 / 3
R_konv_innen3 = 1 / A_i_3 * alpha_i_3
R_konv_aussen3 = 1 / A_a_3 * alpha_a_3
R_waermeleitung3 = np.log(d_aa / d_ai) / (2 * np.pi * l3 * lambda_fluid_3)

R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3
# delta_T2_2 = T3 - T2_sattdampf
cp_fluid_3 = PropsSI('C', 'T', T2_R, 'Q', 0, fluid)
# Q_zu3 = m_ORC * (h3 - h2_sattdampf)
dTA = Thoch_H - Thoch_L

from scipy.optimize import fsolve

T3 = fsolve(solveT3, 350., args=(Q_zu3, d_i, l3, alpha_i_3, T2_sattdampf, dTA))

h3 = PropsSI('H', 'T', T3, 'P', p2, fluid)
Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3


# Turbine

from isentroper_Wirkungsgrad_Expander import isentroper_Wirkungsgrad

n = 5000
eta_Expander = isentroper_Wirkungsgrad(m_ORC, n)
s3 = PropsSI('S', 'P', p2, 'H', h3, fluid)
p4 = p1  # Druckverhältnis variieren
h4 = PropsSI('H', 'S', s3, 'P', p4, fluid)
T4 = PropsSI('T', 'P', p4, 'H', h4, fluid)
w_t = (h4 - h3)
P_t = m_ORC * (h4 - h3) * eta_Expander

# TODO Druckverhältnis implementieren und variieren
verhaeltnis = p2 / p4


# Kondensator 2 SF -> UK

rho_k2 = PropsSI('D', 'T', T4, 'P', p4, fluid)  # kg/m3
v_k2 = m_ORC / rho_k2
c_k2 = v_k2 / A_quer
viscosity_k2 = PropsSI('VISCOSITY', 'T', T4, 'P', p4, fluid)

re_k2 = rho_k2 * c_k2 * d_i / viscosity_k2
lambda_fluid_k2 = PropsSI('CONDUCTIVITY', 'T', T4, 'P', p4, fluid)
pr_k2 = PropsSI('PRANDTL', 'T', T4, 'P', p4, fluid)

alpha_i_k2 = alpha_inside_tube(re_k2, pr_k2, lambda_fluid_k2, d_i)
alpha_a_k2 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k2)

Q_ab2 = m_ORC * (h4_sattdampf - h1)

dTA_k2 = T4_sattdampf - T1
dTB_k2 = 20  # Kühlwasser etc.
l_k2 = Q_ab2 / (np.pi * d_i * alpha_i_k2 * ((dTA_k2) - (dTB_k2) / np.log(dTA_k2 / dTB_k2)))
A_i = 2 * np.pi * d_i / 2 * l_k2 / 3
A_a = 2 * np.pi * d_ai / 2 * l_k2 / 3
R_konv_innen = 1 / A_i * alpha_i_k2
R_konv_aussen = 1 / A_a * alpha_a_k2
R_waermeleitung = np.log(d_aa / d_ai) / (2 * np.pi * l_k2 * lambda_fluid_k2)


"Berechnung thermischer Wirkungsgrad"

P_netto = abs(P_t + P_p)
eta_th = P_netto / (Q_zu_ges)

plt.plot(m_ORC, eta_th, color='black', marker='.', linestyle='-')
plt.xlabel('Massenstrom', fontsize=16)
plt.ylabel('thermischer Wirkungsgrad', fontsize=16)

# plt.legend(loc='best')
plt.show()









