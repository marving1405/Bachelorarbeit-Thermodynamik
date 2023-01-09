# -*- coding: utf-8 -*-
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

fluid = "PROPANE" #REFPROP::
m_ORC = 10E-3 #TODO implemtieren Verhältnis aus ORC Massenstrom und Oelmassenstrom

"""
Pumpe: Zustand 1 so gewählt, dass Fluid bei 1 bar und unterkühlt vorliegt

"""
p1 = 100000 #kPa
T1 = 229 #Kelvin
etaP = 0.9 #wird hier als konstant angesehen

# Berechnung der zu verrichtenden Pumpenarbeit
# Eintrittszustand bei 1 bar unterkühlte Flüssigkeit

h1 = PropsSI('H','T',T1,'P',p1,fluid)
p2 = 2000000 #Pa
v1 = 1 / (PropsSI('D','P',p1,'T',T1,fluid)) #m3/kg

w_p = (v1 * (p2-p1)) #J
P_p = (w_p * m_ORC) / (etaP)
h2 = w_p + h1
T2 = PropsSI('T','H',h2,'P',p2,fluid)


""" Verdampfer
    Betrachtet wird ein Doppelrohrwärmeübertrager
    Als Rohrdurchmesser werden Innen 10mm und außen 14mm festgelegt
    Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
    """

d_i = 10E-3 #m
d_ai = 14E-3
d_aa = 16E-3

T3 = 367.456 #K Berechnet mit LMTD-Methode -> Solver-Funktion #TODO T3 wird erst weiter unten berechnet, später an richtige Stelle
h3 = PropsSI('H','T',T3,'P',p2,fluid) #TODO h3 analog

#TODO Implementieren Massenstromverhältnis
#for x in np.arange(0.1,10,0.1):
m_oel = 4 * m_ORC

cp_oel = 1.9 #kJ/kg*K
Tmittel_H = 100 + 273.15 #K

'''
Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)

'''
T2_siedend = PropsSI('T','P',p2,'Q',0,fluid)
cp_fluid_1 = PropsSI('C', 'T', T2, 'P', p2, fluid)
Q_zu1 = m_ORC * cp_fluid_1 * (T2_siedend - T2)
A_quer = np.pi * (d_i/2)**2 #m2

from if_Abfrage_Rekuperation import abfrage
#abfrage()

rho_1 = PropsSI('D','T',T2,'P',p2,fluid)  #kg/m3
v_1 = m_ORC / rho_1
c_1 = v_1 / A_quer
viscosity_1 = PropsSI('VISCOSITY','T',T2,'P',p2,fluid)


#re_1 = rho_1 * c_1 * d_i / viscosity_1
lambda_fluid_1 = PropsSI('CONDUCTIVITY','T',T2,'P',p2,fluid)
pr_1 = PropsSI('PRANDTL','T',T2,'P',p2,fluid)
#alpha_i_1 = alpha_inside_tube(re_1,pr_1,lambda_fluid_1, d_i)
alpha_a_1 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_1)

from calculate_alpha_aw import alpha_1P_i

alpha_i_1 = alpha_1P_i(p2,T2,fluid,m_ORC,d_i)






Tlow_H = 60 + 273.15 #K
Tlow_L = 30 + 273.15 #K
l1 = Q_zu1 / (np.pi * d_i * alpha_i_1 * (Tlow_H - Tlow_L))

A_i = 2 * np.pi * d_i/2 * l1/3
A_a = 2 * np.pi * d_ai/2 * l1/3
R_konv_innen1 = 1 / A_i * alpha_i_1
R_konv_aussen1 = 1 / A_a * alpha_a_1
R_waermeleitung1 = np.log(d_aa/d_ai) / (2 * np.pi * l1 * lambda_fluid_1)
R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1


'''
Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
isotherme Zustandsänderung, daher über 1.HS
'''

lambda_fluid_2 = PropsSI('CONDUCTIVITY','T',T2_siedend,'Q',0,fluid)
alpha_i_zweiphasig = 600
alpha_a_zweiphasig = 400
# TODO andere alphas aus VDI Waermeatlas ok?


h2_sattdampf = PropsSI('H','P',p2,'Q',1,fluid)
h2_siedend = PropsSI('H','P',p2,'Q',0,fluid)
Q_zu2 = m_ORC * (h2_sattdampf - h2_siedend)
Tmittel_L = Tmittel_H - ((Q_zu2/1000)/(m_oel * cp_oel))

l2 = Q_zu2 / (np.pi * d_i * alpha_i_zweiphasig * (Tmittel_H-Tmittel_L))
A_i_2 = 2 * np.pi * d_i/2 * l2/3
A_a_2 = 2 * np.pi * d_ai/2 * l2/3

R_konv_innen2 = 1 / A_i_2 * alpha_i_zweiphasig
R_konv_aussen2 = 1 / A_a_2 * alpha_a_zweiphasig
R_waermeleitung2 = np.log(d_aa/d_ai) / (2 * np.pi * l2 * lambda_fluid_2)
R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2


'''
Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
'''

#TODO Solver programmieren
Thoch_H = 180 + 273.15 #K
Thoch_L = 110 + 273.15 #K
l3 = 5 #m

h3 = PropsSI('H','T',T3,'P',p2,fluid)
T2_sattdampf = PropsSI('T','P',p2,'Q',1,fluid)

rho_3 = PropsSI('D','T',T2_sattdampf,'Q',1,fluid)  #kg/m3
v_3 = m_ORC / rho_3
c_3 = v_3 / A_quer
viscosity_3 = PropsSI('VISCOSITY','T',T2_sattdampf,'Q',1,fluid)

re_3 = rho_3 * c_3 * d_i / viscosity_3
lambda_fluid_3 = PropsSI('CONDUCTIVITY','T',T2_sattdampf,'Q',1,fluid)
pr_3 = PropsSI('PRANDTL','T',T2_sattdampf,'Q',1,fluid)
alpha_i_3 = alpha_inside_tube(re_3,pr_3,lambda_fluid_3, d_i)
alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_3)

A_i_3 = 2 * np.pi * d_i/2 * l3/3
A_a_3 = 2 * np.pi * d_ai/2 * l3/3
R_konv_innen3 = 1 / A_i_3 * alpha_i_3
R_konv_aussen3 = 1 / A_a_3 * alpha_a_3
R_waermeleitung3 = np.log(d_aa/d_ai) / (2 * np.pi * l3 * lambda_fluid_3)


R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3
delta_T2_2 = T3 - T2_sattdampf
cp_fluid_3 = PropsSI('C', 'T', T2_siedend, 'Q', 0, fluid)
Q_zu3 = m_ORC * (h3 - h2_sattdampf)
dTA = Thoch_H - Thoch_L

'''SOLVER
from scipy.optimize import fsolve
from sympy import *


dTB = symbols('dTB')


sol = solve([eq1, eq2], [dTB])
print(dTB)
'''
'''
from sympy import solve
from sympy.solvers import symbol

from sympy.solvers import nsolve
from sympy.solvers import solveset


dTA = Thoch_H - Thoch_L
dTB = T3 - T2_sattdampf

Q_zu3 / (np.pi * d_i * l3 * alpha_i_3) = (dTA - dTB)/(np.log(dTA/dTB)

dTB = symbol('dTB')
solve(a = (dTA-dTB)/(np.log(dTA/dTB), dTB)
'''
'''
def eq1(dTB):
    return ((dTA - dTB) / (np.log(dTA / dTB)))



def eq2():
    return (Q_zu3 / (np.pi * d_i * l3 * alpha_i_3))

from scipy import optimize
sol = optimize.root(eq1, [0, 0], method='hybr')
sol.dTB

'''


Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3


# Turbine

from isentroper_Wirkungsgrad_Expander import isentroper_Wirkungsgrad
n = 5000
eta_Expander = isentroper_Wirkungsgrad(m_ORC, n)
s3 = PropsSI('S','H',h3,'P',p2,fluid)
p4 = p1 #Druckverhältnis variieren
h4 = PropsSI('H','S',s3,'P',p4,fluid)
T4 = PropsSI('T','P',p4,'H',h4,fluid)
w_t = (h4 - h3) * eta_Expander
P_t = m_ORC * (h4 - h3) * eta_Expander
# TODO Druckverhältnis implementieren und variieren
verhaeltnis = p2/p4


# Kondensator 1

#TODO Länge der Kondensatoren mit LMTD berechnen
l_k1 = 5 #m Annahme
rho_k1 = PropsSI('D','T',T4,'P',p4,fluid)  #kg/m3
v_k1 = m_ORC / rho_k1
c_k1 = v_k1 / A_quer
viscosity_k1 = PropsSI('VISCOSITY','T',T4,'P',p4,fluid)


re_k1 = rho_k1 * c_k1 * d_i / viscosity_k1
lambda_fluid_k1 = PropsSI('CONDUCTIVITY','T',T4,'P',p4,fluid)
pr_k1 = PropsSI('PRANDTL','T',T4,'P',p4,fluid)

alpha_i_k1 = alpha_inside_tube(re_k1,pr_k1,lambda_fluid_k1, d_i)
alpha_a_k1 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k1)

h4_siedend = PropsSI('H','P',p4,'Q',1,fluid)
Q_ab1 = m_ORC * (h4 - h4_siedend)

A_i = 2 * np.pi * d_i/2 * l_k1/2
A_a = 2 * np.pi * d_ai/2 * l_k1
R_konv_innen = 1 / A_i * alpha_i_k1
R_konv_aussen = 1 / A_a * alpha_a_k1
R_waermeleitung = np.log(d_aa/d_ai) / (2 * np.pi * l_k1 * lambda_fluid_k1)



# Kondensator 2
l_k2 = 5 #m
rho_k2 = PropsSI('D','T',T4,'P',p4,fluid)  #kg/m3
v_k2 = m_ORC / rho_k2
c_k2 = v_k2 / A_quer
viscosity_k2 = PropsSI('VISCOSITY','T',T4,'P',p4,fluid)


re_k2 = rho_k2 * c_k2 * d_i / viscosity_k2
lambda_fluid_k2 = PropsSI('CONDUCTIVITY','T',T4,'P',p4,fluid)
pr_k2 = PropsSI('PRANDTL','T',T4,'P',p4,fluid)

alpha_i_k2 = alpha_inside_tube(re_k2,pr_k2,lambda_fluid_k2, d_i)
alpha_a_k2 = alpha_outside_tube(d_ai, d_aa, lambda_fluid_k2)

Q_ab2 = m_ORC * (h4_siedend - h1)

R_konv_innen = 1 / A_i * alpha_i_k2
R_konv_aussen = 1 / A_a * alpha_a_k2
R_waermeleitung = np.log(d_aa/d_ai) / (2 * np.pi * l_k2 * lambda_fluid_k2)




"Berechnung thermischer Wirkungsgrad"

P_netto = abs(P_t + P_p)
eta_th = P_netto / (Q_zu_ges)

plt.plot(m_ORC, eta_th, color='black',marker = '.', linestyle = '-')
plt.xlabel('Massenstrom', fontsize=16)
plt.ylabel('thermischer Wirkungsgrad', fontsize=16)

#plt.legend(loc='best')
plt.show()

'''''
from calculate_alpha_aw import alpha_1P_i

if Rekuperation = true:
    # Rekuperator
    
    # Berechnung nach lmtd-Methode Gleichung nach Temperatur mit l = 5m auflösen (mit Q_zu1)
    T5 = 247.113  # K
else: #Berechnung mit WÜ1 und Tanks
    #Programmierung WÜ

'''''

#TODO Kondensatoren mit Rekuperator in if-Schleife implementieren (und WÜ1)







