"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""
import json, CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import sys
from calculate_alpha_aw import alpha_inside_tube, alpha_outside_tube
from Test_fsolve import solveT3
from calculate_alpha_aw import alpha_1P_i
from calculate_alpha_aw import alpha_boiling
from calculate_alpha_aw import alpha_1P_annulus
from scipy.optimize import fsolve
from Stoffdaten_Oel_Funktionen import lambda_Oel
from Stoffdaten_Oel_Funktionen import cp_Oel

from isentroper_Wirkungsgrad_Expander import isentroper_Wirkungsgrad

CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
#for Thoch_H in np.arange(70 + 273.15,200 + 273.15, 1):
#for p2 in np.arange(100000, 2100000, 10000):
a = []
b = []
c = []
d = []
e = []
f = []
g = []


fluid = "REFPROP::PROPANE" #[0.7]&METHANE[0.3]"
T3 = 400 #K
m_ORC = 40E-3  # kg/s
m_WASSER = 60E-3
cp_WASSER = 4.1819 #kJ/kg*K
m_OEL_2 = 250E-3
m_OEL_3 = 60E-3
m_Kuehlmittel1 = 320E-3
h_g = CP.PropsSI('H', 'P', 101325, 'Q', 1, fluid)
h_liq = CP.PropsSI('H', 'P', 101325, 'Q', 0, fluid)
h_v = h_g - h_liq # Vedampfungsenthalpie



"""
Pumpe: Zustand 1 so anpassen, dass Fluid bei gewünschtem Druck unterkühlt vorliegt
"""
p1 = 100000  # Pa
T1 = 229  # Kelvin
etaP = 0.9  # wird hier als konstant angesehen


h1 = CP.PropsSI('H', 'T', T1, 'P', p1, fluid)
p2 = 2000000  # Pa
v1 = 1 / (CP.PropsSI('D', 'P', p1, 'T', T1, fluid))  # m3/kg

w_p = (v1 * (p2 - p1)) / etaP # J
P_p = (w_p * m_ORC)
h2 = w_p + h1
T2 = CP.PropsSI('T', 'H', h2, 'P', p2, fluid)
if CP.PhaseSI('H', h2, 'P', p2, fluid) == 'twophase':
    raise ValueError("no liquid state at pump outlet!")

""" 
Betrachtet wird ein Doppelrohrwärmeübertrager
Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
"""

d_i = 30E-3  # m
d_ai = 32E-3
d_aa = 38E-3

'''
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)
'''

T2_siedend = CP.PropsSI('T', 'P', p2, 'Q', 0, fluid)
h2_siedend = CP.PropsSI('H', 'P', p2, 'Q', 0, fluid)
Q_zu1 = m_ORC * (h2_siedend - h2)

speicherfluid1 = "REFPROP::WATER"
Tlow_H = T2_siedend + 5 #K
Tlow_L = Tlow_H - ((Q_zu1/1000) / (m_WASSER * cp_WASSER))

p_Tank1 = 100000  # Pa

lambda_Wasser = 0.6
alpha_a_1 = alpha_outside_tube(d_ai,d_aa,lambda_Wasser)
alpha_i_1 = alpha_1P_i(p2, T2, fluid, m_ORC, d_i)

dTA_1 = Tlow_L - T2
dTB_1 = Tlow_H - T2_siedend

lambda_Kupfer = 401
A_i = np.pi * d_i
A_a = np.pi * d_ai
R_konv_innen1 = 1 / (A_i * alpha_i_1)
R_konv_aussen1 = 1 / (A_a * alpha_a_1)
R_waermeleitung1 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
l1 = (Q_zu1 * R_ges1) / ((dTA_1 - dTB_1) / np.log(dTA_1 / dTB_1))

'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''
#arbeitsfluid2 = shell heat transfer oil s2
Tmittel_L = T2_siedend + 5  # 5K pinch point
p_Tank2 = 100000  # Pa
lambda_Oel = 0.129
T2_sattdampf = CP.PropsSI('T', 'P', p2, 'Q', 1, fluid)
viscosity2_liq = CP.PropsSI('VISCOSITY', 'Q', 0, 'P', p2, fluid)
viscosity2_gas = CP.PropsSI('VISCOSITY', 'Q', 1, 'P', p2, fluid)
cp2_liq = CP.PropsSI('C', 'P', p2, 'Q', 0, fluid)
sigma = CP.PropsSI('SURFACE_TENSION', 'P', p2, 'Q', 0, fluid)


h2_sattdampf = CP.PropsSI('H', 'P', p2, 'Q', 1, fluid)
Q_zu2 = m_ORC * (h2_sattdampf - h2_siedend)


rho2_siedend = CP.PropsSI('D', 'Q', 0, 'P', p2, fluid)
rho2_sattdampf = CP.PropsSI('D', 'Q', 1, 'P', p2, fluid)
lambda_fluid_2 = CP.PropsSI('CONDUCTIVITY', 'Q', 0, 'P', p2, fluid)
Te = T2_siedend #Tmittel_H + (((Q_zu2/1000) * np.log(d_ai/d_i)) / (2 * np.pi * 25 * lambda_fluid_2))
dPsat = CP.PropsSI('P', 'T', T2_siedend, 'Q', 0, fluid) - CP.PropsSI('P', 'T', T2_siedend, 'Q', 0, fluid)
alpha_i_zweiphasig = alpha_boiling(m_ORC, 0.5, d_i, rho2_siedend, rho2_sattdampf, viscosity2_liq, viscosity2_gas, lambda_fluid_2, cp2_liq, h_v, sigma, dPsat, Te)
alpha_a_2 = alpha_outside_tube(d_ai, d_aa, lambda_Oel)

cp_Oel = 2.2 # in kJ/kg, konstant angenommen
Tmittel_H = Tmittel_L + ((Q_zu2 / 1000) / (m_OEL_2 * cp_Oel))

dTA_2 = Tmittel_L - T2_siedend
dTB_2 = Tmittel_H - T2_sattdampf

A_i_2 = np.pi * d_i
A_a_2 = np.pi * d_ai

R_konv_innen2 = 1 / (A_i_2 * alpha_i_zweiphasig)
R_konv_aussen2 = 1 / (A_a_2 * alpha_a_2)
R_waermeleitung2 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
l2 = (Q_zu2 * R_ges2) / ((dTA_2 - dTB_2) / np.log(dTA_2 / dTB_2))

'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
'''

Thoch_H = T3 + 5  # K

h3 = CP.PropsSI('H', 'T', T3, 'P', p2, fluid)
Q_zu3 = m_ORC * (h3 - h2_sattdampf)
Thoch_L = Thoch_H - ((Q_zu3 / 1000) / (m_OEL_3 * cp_Oel))

dTA_3 = Thoch_L - T2_sattdampf
dTB_3 = Thoch_H - T3
if dTA_3 < 0:
    raise ValueError("temperature of working fluid is higher then storage fluid")



alpha_i_3 = alpha_1P_i(p2, T2_sattdampf, fluid, m_ORC, d_i)
alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_Oel)


A_i_3 = np.pi * d_i
A_a_3 = np.pi * d_ai
R_konv_innen3 = 1 / (A_i_3 * alpha_i_3)
R_konv_aussen3 = 1 / (A_a_3 * alpha_a_3)
R_waermeleitung3 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3
l3 = (Q_zu3 * R_ges3) / ((dTA_3 - dTB_3) / np.log(dTA_3 / dTB_3))

#T3 = fsolve(solveT3, 450., args=(Q_zu3, R_ges3, Thoch_H, dTA_3))
#print(T3)

Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3

'''
Überprüfung Strömungsgeschwindigkeit
'''
rho3 = CP.PropsSI('D', 'T', T3, 'P', p2, fluid)
v_3 = m_ORC / rho3
c_3 = v_3 / A_i_3
if c_3 > 2:
    raise ValueError("the flow rate of the fluid is too high")

'''
Berechnung Turbine
'''


eta_Expander = isentroper_Wirkungsgrad()
s3 = CP.PropsSI('S', 'P', p2, 'H', h3, fluid)
p4 = p1
h4s = CP.PropsSI('H', 'S', s3, 'P', p4, fluid)
T4s = CP.PropsSI('T', 'P', p4, 'H', h4s, fluid)
w_ts = (h4s - h3)
w_t = w_ts * eta_Expander
P_t = m_ORC * w_t
h4 = eta_Expander * w_ts + h3
T4 = CP.PropsSI("T", "H", h4, "P", p1, fluid)

'''
Kondensator 1, ÜD -> SF, Kühlmedium R23
'''
kuehlmittel1 = "REFPROP::R23"
p_Kuehlmittel1 = 1500000  # Pa


h4_siedend = CP.PropsSI('H', 'P', p4, 'Q', 0, fluid)
T4_siedend = CP.PropsSI('T', 'P', p4, 'H', h4_siedend, fluid)
Q_ab1 = m_ORC * (h4 - h1)
Ta_kuehlmittel1 = T4 - 75 #Temperaturdifferenz anpassen!
ha_kuehlmittel1 = CP.PropsSI('H', 'P', p_Kuehlmittel1, 'T', Ta_kuehlmittel1, kuehlmittel1)
he_kuehlmittel1 = ha_kuehlmittel1 - (Q_ab1 / m_Kuehlmittel1)
Te_kuehlmittel1 = CP.PropsSI('T', 'P', p_Kuehlmittel1, 'H', he_kuehlmittel1, kuehlmittel1)

dTA_k1 = T4 - Ta_kuehlmittel1
dTB_k1 = T1 - Te_kuehlmittel1
lambda_fluid_k1 = CP.PropsSI('CONDUCTIVITY', 'T', Te_kuehlmittel1, 'P', p_Kuehlmittel1, kuehlmittel1)
alpha_i_k1 = alpha_1P_i(p4,T4,fluid,m_ORC,d_i)
alpha_a_k1 = alpha_1P_annulus(p4,Te_kuehlmittel1,kuehlmittel1,m_Kuehlmittel1,d_ai,d_aa)

A_i_k1 = np.pi * d_i
A_a_k1 = np.pi * d_ai
R_konv_innen_k1 = 1 / (A_i_k1 * alpha_i_k1)
R_konv_aussen_k1 = 1 / (A_a_k1 * alpha_a_k1)
R_waermeleitung_k1 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges_k1 = R_konv_innen_k1 + R_konv_aussen_k1 + R_waermeleitung_k1
l_k1 = (Q_ab1 * R_ges_k1) / ((dTA_k1 - dTB_k1) / np.log(dTA_k1 / dTB_k1))

'''
Kondensator 2 SF -> UK, Kühlmedium Methanol

kuehlmittel2 = "REFPROP::METHANOL"
p_Kuehlmittel2 = 100000  # Pa
m_Kuehlmittel2 = 10E-3

Q_ab2 = m_ORC * (h4_siedend - h1)

Ta_kuehlmittel2 = T1 - 5  # pinch point temperature = 20K difference
ha_kuehlmittel2 = CP.PropsSI('H', 'P', p_Kuehlmittel2, 'T', Ta_kuehlmittel2, kuehlmittel2)
he_kuehlmittel2 = ha_kuehlmittel2 - (Q_ab2 / m_Kuehlmittel2)
Te_kuehlmittel2 = CP.PropsSI('T', 'P', p_Kuehlmittel2, 'H', he_kuehlmittel2, kuehlmittel2)

dTA_k2 = T4_siedend - Ta_kuehlmittel2
dTB_k2 = T1 - Te_kuehlmittel2
lambda_fluid_k2 = CP.PropsSI('CONDUCTIVITY', 'T', Te_kuehlmittel2, 'P', p_Kuehlmittel2, kuehlmittel2)
alpha_i_k2 = alpha_1P_i(p4,T4_siedend,fluid,m_ORC,d_i)
alpha_a_k2 = alpha_1P_annulus(p4,Te_kuehlmittel2,kuehlmittel2,m_Kuehlmittel2,d_ai,d_aa)

A_i_k2 = np.pi * d_i
A_a_k2 = np.pi * d_ai
R_konv_innen_k2 = 1 / (A_i_k2 * alpha_i_k2)
R_konv_aussen_k2 = 1 / (A_a_k2 * alpha_a_k2)
R_waermeleitung_k2 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges_k2 = R_konv_innen_k2 + R_konv_aussen_k2 + R_waermeleitung_k2
l_k2 = Q_ab2 / ((1/R_ges_k2) * (dTA_k2 - dTB_k2 / np.log(dTA_k2 / dTB_k2)))
'''
Q_ab_ges = Q_ab1


"Berechnung thermischer Wirkungsgrad"
P_netto = abs(P_t + P_p)
eta_th = P_netto / Q_zu_ges

a.append(eta_th)
b.append(Thoch_H)
c.append(P_netto)
d.append(m_ORC)
e.append(p2)
#f.append(v)

"Entropieberechnung"
Tm1 = (Tlow_H - Tlow_L) / (np.log(Tlow_H / Tlow_L))
Tm2 = (Tmittel_H - Tmittel_L) / (np.log(Tmittel_H / Tmittel_L))
Tm3 = (Thoch_H - Thoch_L) / (np.log(Thoch_H / Thoch_L))
Tmk1 = (Ta_kuehlmittel1 - Te_kuehlmittel1) / (np.log(Ta_kuehlmittel1 / Te_kuehlmittel1))
s_irr = -(Q_zu1 / Tm1 + Q_zu2 / Tm2 + Q_zu3 / Tm3 - Q_ab1 / Tmk1)

g.append(s_irr)



plt.plot(b,c,color='blue')
plt.title("thermischer Wirkungsgrad über TH des Reservoirs", fontsize=14)
plt.xlabel('TH des Reservoirs [K]', fontsize=14)
plt.ylabel('thermischer Wirkungsgrad []', fontsize=14)
plt.grid(True)
plt.show()

print(f"Länge des ersten WÜ: {l1} ")
print(f"Länge des zweiten WÜ: {l2} ")
print(f"Länge des dritten WÜ: {l3} ")
print(f"Länge des Kondensators: {l_k1} ")
print(f"Thoch_H: {Thoch_H} ")
print(f"Thoch_L: {Thoch_L} ")
print(f"Tmittel_H: {Tmittel_H} ")
print(f"Tmittel_L: {Tmittel_L} ")
print(f"Tlow_H: {Tlow_H} ")
print(f"Tlow_L: {Tlow_L} ")


