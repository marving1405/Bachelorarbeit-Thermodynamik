"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""
import json, CoolProp.CoolProp as CP
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

fluid = "REFPROP::PROPANE" #[0.7]&METHANE[0.3]"
# TODO Implementieren Massenstromverhältnis
m_ORC = 10E-3  # kg/s
v = 4 # beschreibt das Verhältnis von Arbeits- zu Prozessfluid
m_OEL = v * m_ORC
cp_oel = 1.9  # kJ/kg*K
h_g = CP.PropsSI('H', 'P', 101325, 'Q', 1, fluid)
h_liq = CP.PropsSI('H', 'P', 101325, 'Q', 0, fluid)
h_v = h_g - h_liq # Vedanpfungsenthalpie



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

h1 = CP.PropsSI('H', 'T', T1, 'P', p1, fluid)
p2 = 2000000  # Pa
v1 = 1 / (CP.PropsSI('D', 'P', p1, 'T', T1, fluid))  # m3/kg

w_p = (v1 * (p2 - p1))  # J
P_p = (w_p * m_ORC) / etaP
h2 = w_p + h1
T2 = CP.PropsSI('T', 'H', h2, 'P', p2, fluid)

""" 
Betrachtet wird ein Doppelrohrwärmeübertrager
Als Rohrdurchmesser werden Innen 10mm und außen 14mm festgelegt
Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
"""

d_i = 10E-3  # m
d_ai = 16E-3
d_aa = 18E-3
A_quer = np.pi * (d_i / 2) ** 2  # m2

"""
Tank 1: Befüllt mit Wasser
"""
arbeitsfluid1 = "REFPROP::WATER"
Tlow_H = 60 + 273.15  # K
Tlow_L = 40 + 273.15  # K
p_Tank1 = 100000  # Pa
m_WASSER = m_OEL
'''
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)
'''
T2_siedend = CP.PropsSI('T', 'P', p2, 'Q', 0, fluid)
h2_siedend = CP.PropsSI('H', 'P', p2, 'Q', 0, fluid)
Q_zu1 = m_ORC * (h2_siedend - h2) #TODO Herstellen einer Verknüpfung zu Tanktemperaturen?

alpha_a_1 = alpha_1P_annulus(p_Tank1, Tlow_L, arbeitsfluid1, m_WASSER, d_ai, d_aa)
alpha_i_1 = alpha_1P_i(p2, T2, fluid, m_ORC, d_i)

dTA_1 = T2_siedend - T2
dTB_1 = Tlow_H - Tlow_L

lambda_Tank1 = CP.PropsSI('CONDUCTIVITY', 'T', Tlow_L, 'P', p_Tank1, arbeitsfluid1)  # Temperatur am Eingang gewählt
A_i = np.pi * d_i
A_a = np.pi * d_ai
R_konv_innen1 = 1 / (A_i * alpha_i_1)
R_konv_aussen1 = 1 / (A_a * alpha_a_1)
R_waermeleitung1 = np.log(d_aa / d_ai) / (2 * np.pi * lambda_Tank1)
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
l1 = Q_zu1 / ((1/R_ges1) * (dTA_1 - dTB_1 / np.log(dTA_1 / dTB_1)))

"""
Tank 2: Befüllt mit Mineralöl
"""
#arbeitsfluid2 = "REFPROP::OIL" #TODO Mineralöl
Tlow_H = 60 + 273.15  # K
Tlow_L = 40 + 273.15  # K
p_Tank2 = 100000  # Pa

'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''
T2_sattdampf = CP.PropsSI('T', 'P', p2, 'Q', 1, fluid)
viscosity2_liq = CP.PropsSI('VISCOSITY', 'Q', 0, 'P', p2, fluid)
viscosity2_gas = CP.PropsSI('VISCOSITY', 'Q', 1, 'P', p2, fluid)
cp2_liq = CP.PropsSI('C', 'T', T2_siedend, 'Q', 0, fluid)
surface_Tension = CP.PropsSI('SURFACE_TENSION', 'P', p2, 'Q', 0, fluid)

#TODO dPsat und Te für alpha boiling berechnen

rho2_siedend = CP.PropsSI('D', 'Q', 0, 'P', p2, fluid) #TODO rhos ändern
rho2_sattdampf = CP.PropsSI('D', 'Q', 1, 'P', p2, fluid)

alpha_i_zweiphasig = 600  # alpha_boiling(m_ORC,0,d_i,rho2_siedend,rho2_sattdampf,viscosity2_liq,viscosity_2_gas,lambda_fluid_2,cp2_liq,h_v,surface_Tension,dPsat,T2_siedend) # TODO alpha-Berechnung zweiphasig
alpha_a_zweiphasig = 400 #TODO Wert annehmen, da Fluid seinen Zustand nicht ändert und nicht in Stoffdatenbank ist (aus VDI-Wärmeatlas)


Tmittel_H = 100 + 273.15  # K
h2_sattdampf = CP.PropsSI('H', 'P', p2, 'Q', 1, fluid)
Q_zu2 = m_ORC * (h2_sattdampf - h2_siedend) #TODO Herstellen einer Verknüpfung zu Tanktemperaturen
Tmittel_L = Tmittel_H - ((Q_zu2 / 1000) / (m_OEL * cp_oel))
lambda_Tank2 = 0.8 #CP.PropsSI('CONDUCTIVITY', 'T', Tmittel_L, 'P', p_Tank2, arbeitsfluid2) #TODO Stoffdaten von Mineralöl

A_i_2 = np.pi * d_i
A_a_2 = np.pi * d_ai

R_konv_innen2 = 1 / (A_i_2 * alpha_i_zweiphasig)
R_konv_aussen2 = 1 / (A_a_2 * alpha_a_zweiphasig)
R_waermeleitung2 = np.log(d_aa / d_ai) / (2 * np.pi * lambda_Tank2)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
l2 = Q_zu2 / ((1/R_ges2) * (Tmittel_H - Tmittel_L))  #Nur Durchschnittstemperatur, da isotherme Zustandsänderung

'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
'''
Thoch_H = 160 + 273.15  # K
Thoch_L = 120 + 273.15  # K
dTA_3 = Thoch_H - Thoch_L
l3 = 20  # m festgelegt
Q_zu3 = m_OEL * cp_oel * (Thoch_H - Thoch_L) * 1000


lambda_fluid_3 = CP.PropsSI('CONDUCTIVITY', 'T', T2_sattdampf, 'Q', 1, fluid)

alpha_i_3 = alpha_1P_i(p2, T2_sattdampf, fluid, m_ORC, d_i)
alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_3) #TODO Stoffdaten Mineralöl

A_i_3 = np.pi * d_i * l3
A_a_3 = np.pi * d_ai * l3
R_konv_innen3 = 1 / (A_i_3 * alpha_i_3)
R_konv_aussen3 = 1 / (A_a_3 * alpha_a_3)
R_waermeleitung3 = np.log(d_aa / d_ai) / (2 * np.pi * l3 * lambda_fluid_3)
R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3
#TODO Gesamtwiderstand in SolveT3 implementieren
T3 = fsolve(solveT3, 350., args=(Q_zu3, R_ges3, T2_sattdampf, dTA_3))
print(R_ges3)
print(T3)
h3 = CP.PropsSI('H', 'T', T3[0], 'P', p2, fluid)
x3 = CP.PropsSI('Q', 'T', T3[0], 'P', p2, fluid)/1000
Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3

'''
Berechnung Turbine
'''

n = 5000
eta_Expander = isentroper_Wirkungsgrad(m_ORC, n)  # TODO isentropen Wirkungsgrad Funktion implementieren
s3 = CP.PropsSI('S', 'P', p2, 'H', h3, fluid)
p4 = p1  # Druckverhältnis variieren
h4 = CP.PropsSI('H', 'S', s3, 'P', p4, fluid)
T4 = CP.PropsSI('T', 'P', p4, 'H', h4, fluid)
w_t = (h4 - h3)
P_t = m_ORC * w_t * eta_Expander

# TODO Druckverhältnis implementieren und variieren
verhaeltnis = p2 / p4

"""
Kondensator 1: Kühlmittel Methanol
"""
kuehlmittel1 = "REFPROP::METHANOL"
p_Kuehlmittel1 = 100000  # Pa
m_Kuehlmittel1 = 40E-3

h4_siedend = CP.PropsSI('H', 'P', p4, 'Q', 0, fluid)
T4_siedend = CP.PropsSI('T', 'P', p4, 'H', h4_siedend, fluid)
Q_ab1 = m_ORC * (h4 - h4_siedend)
Te_kuehlmittel1 = T4 - 30 #pinch point temperature = 30K difference
he_kuehlmittel = CP.PropsSI('H', 'P', p_Kuehlmittel1, 'T', Te_kuehlmittel1, kuehlmittel1)
ha_kuehlmittel = Q_ab1 / (m_Kuehlmittel1) + he_kuehlmittel
Ta_kuehlmittel1 = CP.PropsSI('T', 'P', p_Kuehlmittel1, 'H', ha_kuehlmittel, kuehlmittel1)
print(Te_kuehlmittel1)
print(Ta_kuehlmittel1)
dTA_k1 = T4 - T4_siedend
dTB_k1 = Ta_kuehlmittel1 - Te_kuehlmittel1
print(dTB_k1)


'''
Kondensator 1, ÜD -> SF
'''

lambda_fluid_k1 = CP.PropsSI('CONDUCTIVITY', 'T', T4, 'P', p4, kuehlmittel1)
alpha_i_k1 = alpha_1P_i(p4,T4,fluid,m_ORC,d_i)
alpha_a_k1 = alpha_1P_annulus(p4,T4,kuehlmittel1,m_Kuehlmittel1,d_ai,d_aa)





l_k1 = Q_ab1 / (np.pi * d_i * alpha_i_k1 * ((dTA_k1 - dTB_k1) / (np.log(dTA_k1 / dTB_k1))))
print(l_k1)
A_i_k1 = np.pi * d_i
A_a_k1 = np.pi * d_ai
R_konv_innen_k1 = 1 / (A_i_k1 * alpha_i_k1)
R_konv_aussen_k1 = 1 / (A_a_k1 * alpha_a_k1)
R_waermeleitung_k1 = np.log(d_aa / d_ai) / (2 * np.pi * lambda_fluid_k1)

'''
Kondensator 2 SF -> UK
'''
rho_k2 = CP.PropsSI('D', 'T', T4, 'P', p4, fluid)  # kg/m3
v_k2 = m_ORC / rho_k2
c_k2 = v_k2 / A_quer
viscosity_k2 = CP.PropsSI('VISCOSITY', 'T', T4, 'P', p4, fluid)
re_k2 = rho_k2 * c_k2 * d_i / viscosity_k2
lambda_fluid_k2 = CP.PropsSI('CONDUCTIVITY', 'T', T4, 'P', p4, fluid)
pr_k2 = CP.PropsSI('PRANDTL', 'T', T4, 'P', p4, fluid)
alpha_i_k2 = alpha_inside_tube(re_k2, pr_k2, lambda_fluid_k2, d_i)
alpha_a_k2 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k2)

Q_ab2 = m_ORC * (h4_siedend - h1)
dTA_k2 = T4_siedend - T1
dTB_k2 = 20  # TODO Kühlmittel? und Temperaturdifferenz
l_k2 = Q_ab2 / (np.pi * d_i * alpha_i_k2 * (dTA_k2 - dTB_k2 / np.log(dTA_k2 / dTB_k2)))


A_i_k2 = (2 * np.pi * d_i) / (2 * l_k2)
A_a_k2 = (2 * np.pi * d_ai) / (2 * l_k2)
R_konv_innen_k2 = 1 / (A_i_k2 * alpha_i_k2)
R_konv_aussen_k2 = 1 / (A_a_k2 * alpha_a_k2)
R_waermeleitung_k2 = np.log(d_aa / d_ai) / (2 * np.pi * l_k2 * lambda_fluid_k2)

"Berechnung thermischer Wirkungsgrad"
P_netto = abs(P_t + P_p)
eta_th = P_netto / Q_zu_ges

plt.plot(p2, eta_th, color='black', marker='.', linestyle='-')
plt.xlabel('Verdampfungsdruck', fontsize=16)
plt.ylabel('thermischer Wirkungsgrad', fontsize=16)

# plt.legend(loc='best')
plt.show()
print(eta_th)
