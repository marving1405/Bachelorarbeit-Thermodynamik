"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""
import json, CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from calculate_alpha_aw import alpha_inside_tube, alpha_outside_tube
from Test_fsolve import solveT3
from calculate_alpha_aw import alpha_1P_i
from calculate_alpha_aw import alpha_boiling
from calculate_alpha_aw import alpha_1P_annulus
from scipy.optimize import fsolve
from isentroper_Wirkungsgrad_Expander import isentroper_Wirkungsgrad
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')


fluid = "REFPROP::PROPANE"
# TODO Implementieren Massenstromverhältnis
m_ORC = 10E-3  # kg/s
m_oel = 4 * m_ORC
cp_oel = 1.9  # kJ/kg*K
h_g = PropsSI('H', 'P', 101325, 'Q', 1, fluid)
h_liq = PropsSI('H', 'P', 101325, 'Q', 0, fluid)
h_v = h_g - h_liq

"""
Pumpe: Zustand 1 so gewählt, dass Fluid bei 1 bar und unterkühlt vorliegt
"""
p1 = 100000  # Pa
T1 = 229  # Kelvin
etaP = 0.9  # wird hier als konstant angesehen

"""
Berechnung der zu verrichtenden Pumpenarbeit
Eintrittszustand bei 1 bar unterkühlte Flüssigkeit
"""

h1 = PropsSI('H', 'T', T1, 'P', p1, fluid)
p2 = 2000000  # Pa
v1 = 1 / (PropsSI('D', 'P', p1, 'T', T1, fluid))  # m3/kg

w_p = (v1 * (p2 - p1))  # J
P_p = (w_p * m_ORC) / etaP
h2 = w_p + h1
T2 = PropsSI('T', 'H', h2, 'P', p2, fluid)

""" 
Betrachtet wird ein Doppelrohrwärmeübertrager
Als Rohrdurchmesser werden Innen 10mm und außen 14mm festgelegt
Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
"""

d_i = 10E-3  # m
d_ai = 14E-3
d_aa = 16E-3
A_quer = np.pi * (d_i / 2) ** 2  # m2

'''
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)
'''
T2_siedend = PropsSI('T', 'P', p2, 'Q', 0, fluid)
h2_siedend = PropsSI('H', 'P', p2, 'Q', 0, fluid)
cp_fluid_1 = PropsSI('C', 'T', T2, 'P', p2, fluid)
Q_zu1 = m_ORC * (h2_siedend - h2) #TODO Herstellen einer Verknüpfung zu Tanktemperaturen

rho_1 = PropsSI('D', 'T', T2, 'P', p2, fluid)  # kg/m3
v_1 = m_ORC / rho_1
c_1 = v_1 / A_quer
viscosity_1 = PropsSI('VISCOSITY', 'T', T2, 'P', p2, fluid)
lambda_fluid_1 = PropsSI('CONDUCTIVITY', 'T', T2, 'P', p2, fluid) #TODO außen oder innen?
pr_1 = PropsSI('PRANDTL', 'T', T2, 'P', p2, fluid)
alpha_a_1 = alpha_1P_annulus(p2,T2,fluid,m_oel,d_ai,d_aa)
alpha_i_1 = alpha_1P_i(p2, T2, fluid, m_ORC, d_i)

Tlow_H = 60 + 273.15  # K
Tlow_L = 20 + 273.15  # K
dTA_1 = T2_siedend - T2
dTB_1 = Tlow_H - Tlow_L
l1 = Q_zu1 / (np.pi * d_i * alpha_i_1 * (dTA_1 - dTB_1 / np.log(dTA_1 / dTB_1)))

A_i = 2 * np.pi * d_i / 2 * l1 / 3
A_a = 2 * np.pi * d_ai / 2 * l1 / 3
R_konv_innen1 = 1 / A_i * alpha_i_1
R_konv_aussen1 = 1 / A_a * alpha_a_1
R_waermeleitung1 = np.log(d_aa / d_ai) / (2 * np.pi * l1 * lambda_fluid_1)
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1

'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''
T2_sattdampf = PropsSI('T', 'P', p2, 'Q', 1, fluid)
viscosity2_liq = PropsSI('VISCOSITY', 'Q', 0, 'P', p2, fluid)
viscosity_2_gas = PropsSI('VISCOSITY', 'Q', 1, 'P', p2, fluid)
cp2_liq = PropsSI('C', 'T', T2_siedend, 'Q', 0, fluid)
surface_Tension = PropsSI('SURFACE_TENSION', 'T', T2_siedend, 'Q', 0, fluid)
Te = T2_siedend
p2_sat_Te = PropsSI('P', 'Q', 0, 'T', T2_siedend, fluid)
p2_sat_T = PropsSI('P', 'Q', 1, 'T', T2_sattdampf, fluid)
dPsat = p2_sat_T - p2_sat_Te

rho2_siedend = PropsSI('D', 'T', T2_siedend, 'P', p2, fluid)
rho2_sattdampf = PropsSI('D', 'T', T2_sattdampf, 'P', p2, fluid)
lambda_fluid_2 = PropsSI('CONDUCTIVITY', 'T', T2_siedend, 'Q', 0, fluid)
alpha_i_zweiphasig = 600  # alpha_boiling(m_ORC,0,d_i,rho2_siedend,rho2_sattdampf,viscosity2_liq,viscosity_2_gas,lambda_fluid_2,cp2_liq,h_v,surface_Tension,dPsat,T2_siedend)
alpha_a_zweiphasig = 400
# TODO alpha-Berechnung zweiphasig

Tmittel_H = 100 + 273.15  # K
h2_sattdampf = PropsSI('H', 'P', p2, 'Q', 1, fluid)
Q_zu2 = m_ORC * (h2_sattdampf - h2_siedend) #TODO Herstellen einer Verknüpfung zu Tanktemperaturen
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
Thoch_H = 140 + 273.15  # K
Thoch_L = 120 + 273.15  # K
dTA_3 = Thoch_H - Thoch_L
l3 = 5  # m festgelegt
Q_zu3 = m_oel * cp_oel * (Thoch_H - Thoch_L) * 1000

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

T3 = fsolve(solveT3, 350., args=(Q_zu3, d_i, l3, alpha_i_3, T2_sattdampf, dTA_3))

h3 = PropsSI('H', 'T', T3[0], 'P', p2, fluid)
Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3

'''
Berechnung Turbine
'''

n = 5000
eta_Expander = isentroper_Wirkungsgrad(m_ORC, n)  # TODO isentropen Wirkungsgrad Funktion implementieren
s3 = PropsSI('S', 'P', p2, 'H', h3, fluid)
p4 = p1  # Druckverhältnis variieren
h4 = PropsSI('H', 'S', s3, 'P', p4, fluid)
T4 = PropsSI('T', 'P', p4, 'H', h4, fluid)
w_t = (h4 - h3)
P_t = m_ORC * w_t * eta_Expander

# TODO Druckverhältnis implementieren und variieren
verhaeltnis = p2 / p4

'''
Kondensator 1, ÜD -> SF
'''
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
dTB_k1 = 40  # TODO Kühlmittel? und Temperaturdifferenz
l_k1 = Q_ab1 / (np.pi * d_i * alpha_i_k1 * (dTA_k1 - dTB_k1 / np.log(dTA_k1 / dTB_k1)))

A_i_k1 = 2 * np.pi * d_i / 2 * l_k1 / 2
A_a_k1 = 2 * np.pi * d_ai / 2 * l_k1
R_konv_innen = 1 / A_i * alpha_i_k1
R_konv_aussen = 1 / A_a * alpha_a_k1
R_waermeleitung = np.log(d_aa / d_ai) / (2 * np.pi * l_k1 * lambda_fluid_k1)

'''
Kondensator 2 SF -> UK
'''
rho_k2 = PropsSI('D', 'T', T4, 'P', p4, fluid)  # kg/m3
v_k2 = m_ORC / rho_k2
c_k2 = v_k2 / A_quer
viscosity_k2 = PropsSI('VISCOSITY', 'T', T4, 'P', p4, fluid)
re_k2 = rho_k2 * c_k2 * d_i / viscosity_k2
lambda_fluid_k2 = PropsSI('CONDUCTIVITY', 'T', T4, 'P', p4, fluid)
pr_k2 = PropsSI('PRANDTL', 'T', T4, 'P', p4, fluid)
alpha_i_k2 = alpha_inside_tube(re_k2, pr_k2, lambda_fluid_k2, d_i)
alpha_a_k2 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k2)

Q_ab2 = m_ORC * (h4_siedend - h1)
dTA_k2 = T4_siedend - T1
dTB_k2 = 20  # TODO Kühlmittel? und Temperaturdifferenz
l_k2 = Q_ab2 / (np.pi * d_i * alpha_i_k2 * (dTA_k2 - dTB_k2 / np.log(dTA_k2 / dTB_k2)))

R_konv_innen = 1 / A_i * alpha_i_k2
R_konv_aussen = 1 / A_a * alpha_a_k2
R_waermeleitung = np.log(d_aa / d_ai) / (2 * np.pi * l_k2 * lambda_fluid_k2)

"Berechnung thermischer Wirkungsgrad"
P_netto = abs(P_t + P_p)
eta_th = P_netto / Q_zu_ges

plt.plot(m_ORC, eta_th, color='black', marker='.', linestyle='-')
plt.xlabel('Massenstrom', fontsize=16)
plt.ylabel('thermischer Wirkungsgrad', fontsize=16)

# plt.legend(loc='best')
plt.show()
print(eta_th)
