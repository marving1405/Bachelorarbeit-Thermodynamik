# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""

import numpy as np
from calculate_alpha_aw import alpha_inside_tube, alpha_outside_tube
fluid = "REFPROP::PROPANE"
# Pumpe
"""
Zustand 1 so gewählt, dass Fluid bei 1 bar und unterkühlt vorliegt

"""
p1 = 100000 #kPa
T1 = 229 #Kelvin
etaP = 0.9

# Berechnung der zu verrichtenden Pumpenarbeit
# Eintrittszustand bei 1 bar unterkühlte Flüssigkeit

# Import Stoffdaten
from CoolProp.CoolProp import PropsSI

h1 = PropsSI('H','T',T1,'P',p1,fluid)
p2 = 2000000 #Pa
v1 = 1 / (PropsSI('D','P',p1,'T',T1,fluid)) #m3/kg EINHEITENPROBLEM

#TODO for-Schleife implementieren
# Massenstrom
m = 10E-3 #kg/s
w_p = (v1 * (p2-p1)) #J EINHEITENPROBLEM
P_p = (w_p * m) / (etaP) #W???
h2 = w_p + h1
T2 = PropsSI('T','H',h2,'P',p2,fluid)
        

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
A_quer = np.pi * d_i**2 #m2
rho = PropsSI('D','T',T1,'P',p1,fluid)  #kg/m3
v = m / rho
c = v / A_quer
viscosity = PropsSI('VISCOSITY','T',T2,'P',p2,fluid)


re = rho * c * d_i / viscosity
lambda_fluid = PropsSI('CONDUCTIVITY','T',T2,'P',p2,fluid)
pr = PropsSI('PRANDTL','T',T2,'P',p2,fluid)
alpha_i = alpha_inside_tube(re,pr,lambda_fluid, d_i)


alpha_a = alpha_outside_tube(d_ai, d_aa, lambda_fluid)


# Laenge des Waermeuebertragers
l = 9 #m

#U = Wärmeübertragungskoeffizient
 #m
#A = #Austauschfläche m2
T2_heat = 200 + 273.15
T3_heat = 100 + 273.15
T3 = 423.15 #Mindesttemperatur ist 150°C max. 200°C
#p2 = p3
h3 = PropsSI('H','T',T3,'P',p2,fluid)

delta_Ta = T2_heat - T3_heat
delta_Tb = T3 - T2
#lmtd = (delta_Ta - delta_Tb) / (np.log(delta_Ta/delta_Tb))

'''
# Temperaturen nochmal updaten -> was kühlt wie sehr ab?
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)

'''
A_i = 2 * np.pi * d_i/2 * l/3
A_a = 2 * np.pi * d_ai/2 * l/3
R_konv_innen1 = 1 / A_i * alpha_i
R_konv_aussen1 = 1 / A_a * alpha_a
R_waermeleitung1 = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid)

T2_siedend = PropsSI('T','P',p2,'Q',0,fluid)
delta_T2_1 = T2_siedend - T2
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
Q_zu1 = (1 / R_ges1) * delta_T2_1



'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''

alpha_i_zweiphasig = 600
alpha_a_zweiphasig = 400
# TODO andere alphas aus VDI Waermeatlas ok?
R_konv_innen2 = 1 / A_i * alpha_i_zweiphasig
R_konv_aussen2 = 1 / A_a * alpha_a_zweiphasig
R_waermeleitung2 = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
h2_sattdampf = PropsSI('H','P',p2,'Q',1,fluid)
h2_siedend = PropsSI('H','P',p2,'T',T2_siedend,fluid)
Q_zu2 = h2_sattdampf - h2_siedend

'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)

'''

R_konv_innen3 = 1 / A_i * alpha_i
R_konv_aussen3 = 1 / A_a * alpha_a
R_waermeleitung3 = np.log(d_aa/d_ai) / (2* np.pi * l * lambda_fluid)

h3 = PropsSI('H','T',T3,'P',p2,fluid)
T2_sattdampf = PropsSI('T','P',p2,'Q',1,fluid)

R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3
delta_T2_2 = T3 - T2_sattdampf
Q_zu3 = (1 / R_ges3) * delta_T2_2


Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3


# Turbine

s3 = PropsSI('S','H',h3,'P',p2,fluid)
etaT = 0.8
p4 = p1
h4 = PropsSI('H','S',s3,'P',p4,fluid)
T4 = PropsSI('T','P',p4,'H',h4,fluid)
w_t = (h4 - h3) * etaT
P_t = m * (h4 - h3) * etaT



# Kondensator 1
# Ist hier das alpha gleich wie beim Verdampfer? Und Einheiten der Wärme

A_i = 2 * np.pi * d_i/2 * l
A_a = 2 * np.pi * d_ai/2 * l
R_konv_innen = 1 / A_i * alpha_i
R_konv_aussen = 1 / A_a * alpha_a
R_waermeleitung = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid)

# Kondensator 2
# Ist hier das alpha gleich wie beim Verdampfer? Und Einheiten der Wärme

A_i = 2 * np.pi * d_i/2 * l
A_a = 2 * np.pi * d_ai/2 * l
R_konv_innen = 1 / A_i * alpha_i
R_konv_aussen = 1 / A_a * alpha_a
R_waermeleitung = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid)

delta_T41 = T4 - T1
R_ges = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
Q_ab = (1 / R_ges) * delta_T41


"Berechnung thermischer Wirkungsgrad"

w_netto = abs(w_t + w_p)
eta_th = w_netto / (Q_zu_ges)
print(eta_th)