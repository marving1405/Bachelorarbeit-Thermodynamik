import numpy as np

# Dieses Skript berechnet das unbekannte Temperaturniveau der Sekundärfluide für jeden Wärmeübertrager des Prozesses


def solver_for_WU1(x,Q,R,T,dTx):
    return [Q * R - (((x[0]-T)-(dTx))/(np.log(x[0]-T)-np.log(dTx)))]

def solver_for_WU2(x,Q,R,T,dTx):
    return [Q * R - ((dTx-(x[0]-T))/(np.log(dTx)-np.log(x[0]-T)))]

def solver_for_WU3(x,Q,R,T,dTx):
    return [Q * R - (((x[0]-T)-(dTx))/(np.log(x[0]-T)-np.log(dTx)))]

def solveT_K(Te, Q_ab, R_ges_k, T1, dTA_K):
    return [(Q_ab * R_ges_k) - (((dTA_K)-(T1-Te[0]))/((np.log(dTA_K))-(np.log(T1-Te[0]))))]