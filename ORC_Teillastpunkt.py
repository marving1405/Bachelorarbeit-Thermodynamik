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
from Test_fsolve import solveT
from calculate_alpha_aw import alpha_1P_i
from calculate_alpha_aw import alpha_boiling
from calculate_alpha_aw import alpha_1P_annulus
from scipy.optimize import fsolve
from Stoffdaten_Oel_Funktionen import lambda_Oel
from Stoffdaten_Oel_Funktionen import cp_Oel
from isentroper_Wirkungsgrad_Expander import isentroper_Wirkungsgrad
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
#for p2 in np.arange(500000, 2000000, 10000):
#for Thoch_H in np.arange(160+273.15, 200+273.15, 1):
#for p2 in np.arange(500000, 2100000, 100000):

#for v in np.arange(0.8,5,0.1):
a = []
b = []
c = []
d = []
e = []
f = []
g = []

#for m_ORC in np.arange(5E-3,25E-3,1E-3):
fluid = "REFPROP::PROPANE" #[0.7]&METHANE[0.3]"

m_ORC = 40E-3  # kg/s
v = 1 # beschreibt das Verhältnis von Arbeits- zu Prozessfluid
m_OEL = v * m_ORC
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
Als Rohrdurchmesser werden Innen 10mm und außen 16mm festgelegt
Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
"""

d_i = 10E-3  # m 10,12,18 #TODO anpassen der Durchmesser?
d_ai = 12E-3
d_aa = 18E-3

'''
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)
'''
l1 = 5  # m
T2_siedend = CP.PropsSI('T', 'P', p2, 'Q', 0, fluid)
h2_siedend = CP.PropsSI('H', 'P', p2, 'Q', 0, fluid)
Q_zu1 = m_ORC * (h2_siedend - h2)
speicherfluid1 = "REFPROP::WATER"


#Tlow_L = T2_siedend + 5  # pinch point temperature = 5K
Tlow_H = 80 + 273.15 # K
p_Tank1 = 100000  # Pa
m_WASSER = m_OEL

alpha_a_1 = alpha_1P_annulus(p_Tank1, Tlow_H, speicherfluid1, m_WASSER, d_ai, d_aa)
alpha_i_1 = alpha_1P_i(p2, T2, fluid, m_ORC, d_i)

#dTA_1 = Tlow_L - T2
dTB_1 = Tlow_H - T2_siedend

lambda_Kupfer = 401 #Quelle VDI Wärmeatlas CP.PropsSI('CONDUCTIVITY', 'T', Tlow_H, 'P', p_Tank1, speicherfluid1)  # Temperatur am Eingang gewählt
A_i = np.pi * d_i * l1
A_a = np.pi * d_ai * l1
R_konv_innen1 = 1 / (A_i * alpha_i_1)
R_konv_aussen1 = 1 / (A_a * alpha_a_1)
R_waermeleitung1 = np.log(d_ai / d_i) / (2 * np.pi * l1 * lambda_Kupfer)
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1

#l1 = Q_zu1 / ((1/R_ges1) * (dTA_1 - dTB_1 / np.log(dTA_1 / dTB_1))) #TODO Berechnung optimales dTA und dTB statt l1
Tlow_L = fsolve(solveT, 300., args=(Q_zu1, R_ges1, T2, dTB_1))

'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''
#arbeitsfluid2 = shell heat transfer oil s2
l2 = 60 #m
Tmittel_L = T2_siedend + 5  # K
p_Tank2 = 100000  # Pa
lambda_oel_Tmittel_L = lambda_Oel(Tmittel_L)
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
alpha_a_2 = alpha_outside_tube(d_ai, d_aa, lambda_oel_Tmittel_L)

cp_oel_Tmittel_L = cp_Oel(Tmittel_L) # in kJ/kg
#Tmittel_H = Tmittel_L + ((Q_zu2 / 1000) / (m_OEL * cp_oel_Tmittel_L))

A_i_2 = np.pi * d_i * l2
A_a_2 = np.pi * d_ai * l2

dTA_2 = Tmittel_L - T2_siedend
#dTB_2 = Tmittel_H - T2_sattdampf

R_konv_innen2 = 1 / (A_i_2 * alpha_i_zweiphasig)
R_konv_aussen2 = 1 / (A_a_2 * alpha_a_2)
R_waermeleitung2 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
Tmittel_H = fsolve(solveT, 350., args=(Q_zu2, R_ges2, T2_sattdampf, dTA_2))
#l2 = Q_zu2 / ((1/R_ges2) * (dTA_2 - dTB_2) / (np.log(dTA_2 / dTB_2)))

'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
'''
l3 = 90  # 30m festgelegt
Thoch_H = 200 + 273.15  # K
Thoch_L = T2_sattdampf + 5  # K

dTA_3 = Thoch_L - T2_sattdampf
if dTA_3 < 0:
    raise ValueError("temperature of working fluid is higher then storage fluid")
cp_oel_Thoch_H = cp_Oel(Thoch_H)
Q_zu3 = m_OEL * cp_oel_Thoch_H * (Thoch_H - Thoch_L) * 1000 #W

lambda_oel_Thoch_H = lambda_Oel(Thoch_H)
alpha_i_3 = alpha_1P_i(p2, T2_sattdampf, fluid, m_ORC, d_i)
alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_oel_Thoch_H)

A_i_3 = np.pi * d_i * l3
A_a_3 = np.pi * d_ai * l3
R_konv_innen3 = 1 / (A_i_3 * alpha_i_3)
R_konv_aussen3 = 1 / (A_a_3 * alpha_a_3)
R_waermeleitung3 = np.log(d_ai / d_i) / (2 * np.pi * l3 * lambda_Kupfer)
R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3

T3 = fsolve(solveT3, 400., args=(Q_zu3, R_ges3, Thoch_H, dTA_3))
if T3 < T2_sattdampf:
    #break
    raise ValueError('T3 is lower then T2_sattdampf')

h3 = CP.PropsSI('H', 'T', T3[0], 'P', p2, fluid)


'''
Überprüfung Strömungsgeschwindigkeit
'''
rho3 = CP.PropsSI('D', 'T', T3, 'P', p2, fluid)
v_3 = m_ORC / rho3
c_3 = v_3 / A_i_3
if c_3 > 2:
    raise ValueError("the flow rate of the fluid is too high")



Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3

'''
Berechnung Expander
'''

n = 5000
eta_Expander = isentroper_Wirkungsgrad(m_ORC, n)
s3 = CP.PropsSI('S', 'P', p2, 'H', h3, fluid)
p4 = p1
h4s = CP.PropsSI('H', 'S', s3, 'P', p4, fluid)
T4s = CP.PropsSI('T', 'P', p4, 'H', h4s, fluid)
w_ts = (h4s - h3)
w_t = w_ts * eta_Expander
P_t = m_ORC * w_t
h4 = eta_Expander * w_ts + h3 #TODO h4 musste neu berechnet werden
T4 = CP.PropsSI("T", "H", h4, "P", p1, fluid)

'''
Kondensator 1, ÜD -> SF, Kühlmedium Methanol
'''
kuehlmittel1 = "REFPROP::METHANOL"
p_Kuehlmittel1 = 100000  # Pa
m_Kuehlmittel1 = 200E-3

h4_siedend = CP.PropsSI('H', 'P', p4, 'Q', 0, fluid)
T4_siedend = CP.PropsSI('T', 'P', p4, 'H', h4_siedend, fluid)
Q_ab1 = m_ORC * (h4 - h4_siedend)
Ta_kuehlmittel1 = T4_siedend - 5 #pinch point temperature = 20K difference
ha_kuehlmittel1 = CP.PropsSI('H', 'P', p_Kuehlmittel1, 'T', Ta_kuehlmittel1, kuehlmittel1)
he_kuehlmittel1 = ha_kuehlmittel1 - (Q_ab1 / m_Kuehlmittel1)
Te_kuehlmittel1 = CP.PropsSI('T', 'P', p_Kuehlmittel1, 'H', he_kuehlmittel1, kuehlmittel1)

dTA_k1 = T4 - Ta_kuehlmittel1
dTB_k1 = T4_siedend - Te_kuehlmittel1
lambda_fluid_k1 = CP.PropsSI('CONDUCTIVITY', 'T', Te_kuehlmittel1, 'P', p_Kuehlmittel1, kuehlmittel1)
alpha_i_k1 = alpha_1P_i(p4,T4,fluid,m_ORC,d_i)
alpha_a_k1 = alpha_1P_annulus(p4,Te_kuehlmittel1,kuehlmittel1,m_Kuehlmittel1,d_ai,d_aa)

A_i_k1 = np.pi * d_i
A_a_k1 = np.pi * d_ai
R_konv_innen_k1 = 1 / (A_i_k1 * alpha_i_k1)
R_konv_aussen_k1 = 1 / (A_a_k1 * alpha_a_k1)
R_waermeleitung_k1 = np.log(d_ai / d_i) / (2 * np.pi * lambda_Kupfer)
R_ges_k1 = R_konv_innen_k1 + R_konv_aussen_k1 + R_waermeleitung_k1
l_k1 = Q_ab1 / ((1/R_ges_k1) * ((dTA_k1 - dTB_k1) / (np.log(dTA_k1 / dTB_k1))))

'''
Kondensator 2 SF -> UK, Kühlmedium Methanol
'''
kuehlmittel2 = "REFPROP::METHANOL"
p_Kuehlmittel2 = 100000  # Pa
m_Kuehlmittel2 = 20E-3

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
R_waermeleitung_k2 = np.log(d_aa / d_ai) / (2 * np.pi * lambda_Kupfer)
R_ges_k2 = R_konv_innen_k2 + R_konv_aussen_k2 + R_waermeleitung_k2
l_k2 = Q_ab2 / ((1/R_ges_k2) * (dTA_k2 - dTB_k2 / np.log(dTA_k2 / dTB_k2)))

Q_ab_ges = Q_ab1 + Q_ab2



"Berechnung thermischer Wirkungsgrad"
P_netto = abs(P_t + P_p)
eta_th = P_netto / Q_zu_ges
a.append(eta_th)
b.append(Thoch_H)
c.append(P_netto)
d.append(m_ORC)

"Entropieberechnung"
Tm1 = (Tlow_H - Tlow_L) / (np.log(Tlow_H / Tlow_L))
Tm2 = (Tmittel_H - Tmittel_L) / (np.log(Tmittel_H / Tmittel_L))
Tm3 = (Thoch_H - Thoch_L) / (np.log(Thoch_H / Thoch_L))
Tmk1 = (Ta_kuehlmittel1 - Te_kuehlmittel1) / (np.log(Ta_kuehlmittel1 / Te_kuehlmittel1))
Tmk2 = (Ta_kuehlmittel2 - Te_kuehlmittel2) / (np.log(Ta_kuehlmittel2 / Te_kuehlmittel2))
s_irr = -(Q_zu1 / Tm1 + Q_zu2 / Tm2 + Q_zu3 / Tm3 - Q_ab1 / Tmk1 - Q_ab2 / Tmk2)

e.append(s_irr)
f.append(p2)
g.append(v)

"Überprüfung der Einhaltung des 1. HS"
if Q_ab_ges + P_netto > Q_zu_ges:
    raise SystemExit

    #plt.plot(Thoch_H,eta_th, marker='.',color='black')
    #plt.plot(Thoch_H, P_netto, marker='.', color='black')
    #plt.scatter(p2, P_netto)
    #plt.legend(loc='best')
    #plt.plot(p2, eta_th, color='black',marker='.', linestyle='solid')
    plt.plot(d,c,color='blue')
    plt.title("thermischer Wirkungsgrad über TH des Reservoirs", fontsize=14)
    plt.xlabel('TH des Reservoirs [K]', fontsize=14)
    plt.ylabel('thermischer Wirkungsgrad []', fontsize=14)
    #plt.title("Nettoleistung über TH des Reservoirs", fontsize=14)
    #plt.xlabel('TH des Reservoirs [K]', fontsize=14)
    #plt.ylabel('Nettoleistung [W]', fontsize=14)
    plt.grid(True)
    plt.show()




