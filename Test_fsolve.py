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

def solveT3(T3):

    return [Q_zu3 / (np.pi * d_i * l3 * alpha_i_3) - (dTA - (T3[0] - T2_sattdampf))/(np.log(dTA/(T3[0] - T2_sattdampf)))]

T3 = (fsolve(solveT3, 500.))

print(T3)

# Q_zu3 / (np.pi * d_i * l3 * alpha_i_3) - (dTA - (T3[0] - T2_sattdampf))/(np.log(dTA/(T3[0] - T2_sattdampf))) = 0
