from scipy.optimize import fsolve
import numpy as np
'''''
Q_zu3 = 5320 #W
d_i = 0.01
l3 = 5
alpha_i_3 = 727
dTA = 70
T2_sattdampf = 330
'''''
def solveT3(T3, Q_zu3, R_ges3, T2_sattdampf, dTA):

    return [Q_zu3 / (1/R_ges3) - (dTA - (T3[0] - T2_sattdampf))/(np.log(dTA/(T3[0] - T2_sattdampf)))]

# Q_zu3 / (np.pi * d_i * l3 * alpha_i_3) - (dTA - (T3[0] - T2_sattdampf))/(np.log(dTA/(T3[0] - T2_sattdampf))) = 0

def solveT2_Rekuperator(T2_Rekuperator, Q_abR, d_i, lR, alpha_R, T2, dTAR):

    return [Q_abR / (np.pi * d_i * lR * alpha_R) - (dTAR - (T2_Rekuperator[0] - T2))/(np.log(dTAR/(T2_Rekuperator[0] - T2)))]

