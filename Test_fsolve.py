from scipy.optimize import fsolve
import numpy as np

def solveT3(T3, Q_zu3, R_ges3, Thoch_H, dTA_3):

    return [(Q_zu3 * R_ges3) - ((dTA_3 - (Thoch_H - T3[0]))/(np.log(dTA_3) - np.log(Thoch_H - T3[0])))]

# Q_zu3 / (np.pi * d_i * l3 * alpha_i_3) - (dTA - (T3[0] - T2_sattdampf))/(np.log(dTA/(T3[0] - T2_sattdampf))) = 0

def solveT2_Rekuperator(T2_Rekuperator, Q_abR, d_i, lR, alpha_R, T2, dTAR):

    return [Q_abR / (np.pi * d_i * lR * alpha_R) - (dTAR - (T2_Rekuperator[0] - T2))/(np.log(dTAR/(T2_Rekuperator[0] - T2)))]

def solveT(Tlow_L1 ,Q_zu1, R_ges1, T2, dTB_1):

    return [(Q_zu1 * R_ges1) - ((Tlow_L1[0] - T2) - dTB_1) / np.log((Tlow_L1[0] - T2) / dTB_1)]

def solveT_K(Te_kuehlmittel ,Q_ab1, R_ges_k, T4_siedend, dTB_K):

    return [(Q_ab1 * R_ges_k) - (((T4_siedend-Te_kuehlmittel[0])-dTB_K) / np.log((T4_siedend-Te_kuehlmittel[0]) / (dTB_K)))]