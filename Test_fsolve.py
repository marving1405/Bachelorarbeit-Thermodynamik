from scipy.optimize import fsolve
import numpy as np
import math
Q_zu3 = 5320 #W
d_i = 0.01
l3 = 5
alpha_i_3 = 654

dTA = 70
#dTB = T3 - T2_sattdampf
T2_sattdampf = 330

def equations(x):
    return [51.7862 - ((400 - x)/(math.log(70) - math.log(x - 330+100000)))]

root = fsolve(equations, [1])

#print(T3)

