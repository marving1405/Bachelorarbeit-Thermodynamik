from scipy.optimize import fsolve
import numpy as np

def solveT3(T3, Q_zu3, R_ges3, Thoch_H, dTA_3):

    return [(Q_zu3 * R_ges3) - ((dTA_3 - (Thoch_H - T3[0]))/(np.log(dTA_3) - np.log(Thoch_H - T3[0])))]

# Q_zu3 / (np.pi * d_i * l3 * alpha_i_3) - (dTA - (T3[0] - T2_sattdampf))/(np.log(dTA/(T3[0] - T2_sattdampf))) = 0



def solver_for_WU1(x,Q,R,T,dTx):
    return [Q * R - (((x[0]-T)-(dTx))/(np.log(x[0]-T)-np.log(dTx)))]


def solver_for_WU2(x,Q,R,T,dTx):
    return [Q * R - ((dTx-(x[0]-T))/(np.log(dTx)-np.log(x[0]-T)))]

def solver_for_WU3(x,Q,R,T,dTx):
    return [Q * R - (((x[0]-T)-(dTx))/(np.log(x[0]-T)-np.log(dTx)))]

def solveT(T ,Q_zu, R_ges, T_ref, dT):

    return [(Q_zu * R_ges) - (((T[0] - T_ref) - dT) / (np.log((T[0] - T_ref) - np.log(dT))))]




def solveT_K(Te, Q_ab, R_ges_k, T1, dTA_K):

    return [(Q_ab * R_ges_k) - (((dTA_K)-(T1-Te[0]))/((np.log(dTA_K))-(np.log(T1-Te[0]))))]