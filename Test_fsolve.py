from scipy.optimize import fsolve
import numpy as np
import math
Q_zu3 = 904
d_i = 0.01
l3 = 5
alpha_i_3 = 654

dTA = 70
#dTB = T3 - T2_sattdampf
T2_sattdampf = 330

def equations(vars):
    T3 = vars
    eq1 = Q_zu3 / (np.pi * d_i * l3 * alpha_i_3)
    eq2 = (dTA - (T3 - T2_sattdampf))/(math.log(dTA/(T3 - T2_sattdampf)))
    return [eq1, eq2]


T3 = fsolve(equations, 1)

print(T3)

