# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""

import numpy as np
from calculate_alpha_aw import alpha_inside_tube, alpha_outside_tube

# Pumpe
"""
Zustand 1 so gewählt, dass Fluid bei 1 bar und unterkühlt vorliegt

"""
p1 = 100 #kPa
T1 = 230.74 #Kelvin
etaP = 0.9

# Berechnung der zu verrichtenden Pumpenarbeit
# Eintrittszustand bei 1 bar unterkühlte Flüssigkeit

# Import Stoffdaten
from CoolProp.CoolProp import PropsSI

h1 = PropsSI('H','T',T1,'P',p1,'Propane')
p2 = 2000 #kPa
v1 = 1 / 581.23 #PropsSI('V','P',p1,'T',T1,'Propane')

# Massenstrom
m = 10E-3 #kg/s
w_p = ((v1 * (p2-p1)) / (etaP))
P_p = (w_p * m) #kW
h2 = w_p + h1
T2 = PropsSI('T','H',h2,'P',p2,'Propane')
        

# Verdampfer WÜ

"""
    delta_Ta ist die Temperaturdifferenz des heißen Fluides (Mineralöl)
    delta_Tb ist die Temperaturdifferenz des ORC-Prozessfluides
    Betrachtet wird ein Doppelrohrwärmeübertrager
    Als Rohrdurchmesser werden Innen 10mm und außen 14mm festgelegt
    Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
    """
d_i = 10E-3 #m
d_ai = 14E-3
d_aa = 16E-3

'''
Parameter für alpha-Berechnung
'''
v = 10 #m/s
dyn_viskositaet = 100 # zw. 80-120
lambda_oel = 0.1 #W/m*K
cp_oel = 1.9 #J/g*K
rho = 880 #kg/m3
re = rho * v * d_i / dyn_viskositaet
pr = (dyn_viskositaet * cp_oel) / (lambda_oel)
pr = 2000 # noch zu berechnen
alpha_i = alpha_inside_tube(re,pr,lambda_oel, d_i)

lambda_fluid = 0.03
alpha_a = alpha_outside_tube(d_ai, d_aa, lambda_fluid)

#U = Wärmeübertragungskoeffizient
l = 10 #m
#A = #Austauschfläche m2
T2_heat = 200 + 273.15
T3_heat = 100 + 273.15
T3 = 423.15 #Mindesttemperatur ist 150°C

delta_Ta = T2_heat - T3_heat
delta_Tb = T3 - T2
lmtd = (delta_Ta - delta_Tb) / (np.log(delta_Ta/delta_Tb))

'''
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu Sattdampf)

'''
A_i = 2 * np.pi * d_i/2 * l
A_a = 2 * np.pi * d_ai/2 * l
R_konv_innen1 = 1 / A_i * alpha_i
R_konv_aussen1 = 1 / A_a * alpha_a
R_waermeleitung1 = np.log(d_aa/d_ai) / (2* np.pi * l * lambda_fluid)
T_sattdampf = PropsSI('T','P',p2,'Q',1,'Propane')
delta_T1 = T2 - T_sattdampf
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
Q_zu1 = (1 / R_ges1) * delta_T1
'''
Auslegung des Wärmeübertragers 2 (Sattdampf zu überhitzten Dampf)

'''
h3_heat = PropsSI('H','T',T3_heat,'P',p2,'Propane')

R_konv_innen2 = 1 / A_i * alpha_i
R_konv_aussen2 = 1 / A_a * alpha_a
R_waermeleitung2 = np.log(d_aa/d_ai) / (2* np.pi * l * lambda_fluid)
T_uedampf = PropsSI('T','H',h3_heat,'P',p2,'Propane')
delta_T2 = T_uedampf - T_sattdampf
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
Q_zu2 = (1 / R_ges2) * delta_T2

#Q_zu = U * A * lmtd


# Turbine

etaT = 0.9
p2 = p3
p4 = p1
#T4 =
P_t = m * (h4 - h3) * etaT

# Berechnung von h4 im isentropen Fall



# Kondensator WÜ