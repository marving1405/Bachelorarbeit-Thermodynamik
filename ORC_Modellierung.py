# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""
import json, CoolProp.CoolProp as CP
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')
from CoolProp.CoolProp import PropsSI


# m = list(range(5E-3,40E-3,1E-3)) #kg/s
import numpy as np
import matplotlib.pyplot as plt


from calculate_alpha_aw import alpha_inside_tube, alpha_outside_tube
#for T3 in np.arange(423.15,473.15,1):
fluid = "REFPROP::PROPANE"
#fluid = "PROPANE"


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
#import CoolProp.CoolProp as CP
#from CoolProp.CoolProp import PropsSI

h1 = PropsSI('H','T',T1,'P',p1,fluid)
p2 = 2000000 #Pa
v1 = 1 / (PropsSI('D','P',p1,'T',T1,fluid)) #m3/kg EINHEITENPROBLEM


# Massenstrom

m = 10E-3
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


# Laenge des Waermeuebertragers
l = 10 #m

#U = Wärmeübertragungskoeffizient
 #m
#A = #Austauschfläche m2
Te_1 = 200 + 273.15
#Ta_1 = 100 + 273.15
T3 = 423.15 #Mindesttemperatur ist 150°C max. 200°C
#p2 = p3
h3 = PropsSI('H','T',T3,'P',p2,fluid)

"""
Berechnung der Speicherfluideigenschaften für den Mitteltemperaturtank
"""
m_oel = 40E-3 #kg/s
cp_oel = 1.9 #kJ/kg*K
Tmittel_H = 100 + 273.15 #K

'''
# Temperaturen nochmal updaten -> was kühlt wie sehr ab?
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)

'''

'''
Parameter für alpha-Berechnung für Wärmeübertrager 1 (Rekuperator)
'''
A_quer = np.pi * (d_i/2)**2 #m2
rho_1 = PropsSI('D','T',T2,'P',p2,fluid)  #kg/m3
v_1 = m / rho_1
c_1 = v_1 / A_quer
viscosity_1 = PropsSI('VISCOSITY','T',T2,'P',p2,fluid)


re_1 = rho_1 * c_1 * d_i / viscosity_1
lambda_fluid_1 = PropsSI('CONDUCTIVITY','T',T2,'P',p2,fluid)
pr_1 = PropsSI('PRANDTL','T',T2,'P',p2,fluid)
alpha_i_1 = alpha_inside_tube(re_1,pr_1,lambda_fluid_1, d_i)
alpha_a_1 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_1)

A_i = 2 * np.pi * d_i/2 * l/3
A_a = 2 * np.pi * d_ai/2 * l/3
R_konv_innen1 = 1 / A_i * alpha_i_1
R_konv_aussen1 = 1 / A_a * alpha_a_1
R_waermeleitung1 = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid_1)

T2_siedend = PropsSI('T','P',p2,'Q',0,fluid)
delta_T2_1 = T2_siedend - T2
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
cp_fluid_1 = PropsSI('CPMASS', 'T', T2, 'P', p2, fluid)
Q_zu1 = m * cp_fluid_1 * delta_T2_1


lmtd1 = Q_zu1 / alpha_a_1 * A_a
from sympy.solvers import solve
from sympy import Symbol

#Ta_1 = Symbol('Ta_1')
#solve(lmtd1 = (Te_1 - Ta_1 - delta_T2_1) / (np.log(Te_1 - Ta_1 / delta_T2_1)), Ta_1)
#Q_zu1 = Q_ab1
#delta_T1 = Te_1 - Ta_1
#lmtd1 = (delta_T1 - delta_T2_1) / (np.log(delta_T1 / delta_T2_1))

'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''
# TODO keine Temperaturdifferenz -> inwiefern therm. Widerstände einbinden?
lambda_fluid_2 = PropsSI('CONDUCTIVITY','T',T2_siedend,'P',p2,fluid)
alpha_i_zweiphasig = 600
alpha_a_zweiphasig = 400
# TODO andere alphas aus VDI Waermeatlas ok?
R_konv_innen2 = 1 / A_i * alpha_i_zweiphasig
R_konv_aussen2 = 1 / A_a * alpha_a_zweiphasig
R_waermeleitung2 = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid_2)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
h2_sattdampf = PropsSI('H','P',p2,'Q',1,fluid)/1000
h2_siedend = PropsSI('H','P',p2,'T',T2_siedend,fluid)/1000
Q_zu2 = m * (h2_sattdampf - h2_siedend)
Tmittel_L = Tmittel_H - (Q_zu2/(m_oel * cp_oel))


'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
'''
'''
Parameter für alpha-Berechnung für Wärmeübertrager 3
'''
Thoch_H = 180 + 273.15 #K
Thoch_L = 110 + 273.15 #K


h3 = PropsSI('H','T',T3,'P',p2,fluid)
T2_sattdampf = PropsSI('T','P',p2,'Q',1,fluid)

rho_3 = PropsSI('D','T',T2_sattdampf,'P',p2,fluid)  #kg/m3
v_3 = m / rho_3
c_3 = v_3 / A_quer
viscosity_3 = PropsSI('VISCOSITY','T',T2_sattdampf,'P',p2,fluid)


re_3 = rho_3 * c_3 * d_i / viscosity_3
lambda_fluid_3 = PropsSI('CONDUCTIVITY','T',T2_sattdampf,'P',p2,fluid)
pr_3 = PropsSI('PRANDTL','T',T2_sattdampf,'P',p2,fluid)
alpha_i_3 = alpha_inside_tube(re_3,pr_3,lambda_fluid_3, d_i)
alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_3)


R_konv_innen3 = 1 / A_i * alpha_i_3
R_konv_aussen3 = 1 / A_a * alpha_a_3
R_waermeleitung3 = np.log(d_aa/d_ai) / (2* np.pi * l * lambda_fluid_3)


R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3
delta_T2_2 = T3 - T2_sattdampf
cp_fluid_3 = PropsSI('CPMASS', 'T', T2_siedend, 'P', p2, fluid)
Q_zu3 = m * cp_fluid_3 * delta_T2_2





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
rho_k1 = PropsSI('D','T',T4,'P',p4,fluid)  #kg/m3
v_k1 = m / rho_k1
c_k1 = v_k1 / A_quer
viscosity_k1 = PropsSI('VISCOSITY','T',T4,'P',p4,fluid)


re_k1 = rho_k1 * c_k1 * d_i / viscosity_k1
lambda_fluid_k1 = PropsSI('CONDUCTIVITY','T',T4,'P',p4,fluid)
pr_k1 = PropsSI('PRANDTL','T',T4,'P',p4,fluid)

alpha_i_k1 = alpha_inside_tube(re_k1,pr_k1,lambda_fluid_k1, d_i)
alpha_a_k1 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k1)

A_i = 2 * np.pi * d_i/2 * l/2
A_a = 2 * np.pi * d_ai/2 * l
R_konv_innen = 1 / A_i * alpha_i_k1
R_konv_aussen = 1 / A_a * alpha_a_k1
R_waermeleitung = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid_k1)




# Kondensator 2
# Ist hier das alpha gleich wie beim Verdampfer? Und Einheiten der Wärme
rho_k2 = PropsSI('D','T',T4,'P',p4,fluid)  #kg/m3
v_k2 = m / rho_k2
c_k2 = v_k2 / A_quer
viscosity_k2 = PropsSI('VISCOSITY','T',T4,'P',p4,fluid)


re_k2 = rho_k2 * c_k2 * d_i / viscosity_k2
lambda_fluid_k2 = PropsSI('CONDUCTIVITY','T',T4,'P',p4,fluid)
pr_k2 = PropsSI('PRANDTL','T',T4,'P',p4,fluid)

alpha_i_k2 = alpha_inside_tube(re_k2,pr_k2,lambda_fluid_k2, d_i)
alpha_a_k2 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k2)

R_konv_innen = 1 / A_i * alpha_i_k2
R_konv_aussen = 1 / A_a * alpha_a_k2
R_waermeleitung = np.log(d_aa/d_ai) / (2 * np.pi * l * lambda_fluid_k2)

delta_T41 = T4 - T1
R_ges = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
cp_fluid_41 = PropsSI('CPMASS', 'T', T4, 'P', p1, fluid)
Q_ab = m * cp_fluid_41 * delta_T41




"Berechnung thermischer Wirkungsgrad"

P_netto = abs(P_t + P_p)
eta_th = P_netto / (Q_zu_ges)
#print(eta_th)

#m = list(range(5E-3,40E-3,1E-3)) #kg/s
plt.plot(T3, eta_th, color='black',marker = '.', linestyle = '-')
plt.xlabel('T3', fontsize=16)
plt.ylabel('eta_th', fontsize=16)

#plt.legend(loc='best')
plt.show()





