# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:24:32 2022

@author: marvi
"""

import CoolProp.CoolProp as CP

# Pumpe
p1 = 100000
#T1 = input('Bitte geben Sie die Starttemperatur ein: ')
etaP = float(input('Bitte geben Sie den isentropen Pumpenwirkungsgrad ein: '))

# Berechnung der zu verrichtenden Verdichterarbeit
# Eintrittszustand bei 1 bar überhitzter Dampf

# Import Stoffdaten
from CoolProp.CoolProp import PropsSI
T1_min = -42.38 + 273.15
#T1 = PropsSI('T','P',p1,'Q',0,'Propane') unbekannt
h1 = PropsSI('H','T',T1_min,'P',p1,'Propane')
p2 = int(input('Bitte geben Sie den Enddruck in Zustand 2 bis maximal 20bar an: '))
T3 = int(input('Bitte geben Sie die Temperatur in Zustand 3 zwischen 150-200°C an: '))
h2 = PropsSI('H','P',p2,'T',T2,'Propane')
# Massenstrom
m = float(input('Bitte geben Sie den Massenstrom zwischen 5-40 g/s ein: '))
wp = ((h2 - h1) / (etaP))
P = (wp * m)

print(P)

#zustand = True
#if zustand = false:
        

# Verdampfer WÜ




# Turbine

# Kondensator WÜ