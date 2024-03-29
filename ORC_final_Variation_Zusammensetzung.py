"""
Created on Thu Nov  3 16:24:32 2022
@author: Marvin Gertenbach
Dieses Skript untersucht, inwiefern sich der Einfluss eines binären zeotropen Fluidgemisches auf die Wärmeübertragung
auswirkt und wie sich dies auf den thermischen Wirkungsgrad und die erzeugte Entropie auswirkt
"""
import json, CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
from calculate_alpha_aw import alpha_outside_tube
from Solver_fsolve import solveT_K
from Solver_fsolve import solver_for_WU1
from Solver_fsolve import solver_for_WU2
from Solver_fsolve import solver_for_WU3
from calculate_alpha_aw import alpha_1P_i
from calculate_alpha_aw import alpha_boiling
from calculate_alpha_aw import alpha_1P_annulus
from scipy.optimize import fsolve
from Stoffdaten_Oel_Funktionen import lambda_Oel
from Stoffdaten_Oel_Funktionen import cp_Oel
from isentroper_Wirkungsgrad_Expander import isentroper_Wirkungsgrad
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')

a = []
b = []
c = []
d = []
e = []
f = []
g = []
h = []
m = []
j = []
k = []
l = []
s = []
t = []
u = []
x_a = []

T_pinch = 5 # K
for x_a in np.arange(0, 1.1, 0.1):
    fluid = f"REFPROP::PROPANE[{round(1-x_a,1)}]&ISOBUTANE[{round(x_a,1)}]"

    l1 = 29 # m
    l2 = 87.5  # m
    l3 = 68 # m
    l_k1 = 22 #69.6  # m
    m_Kuehlmittel1 = 88E-3

    m_ORC = 40E-3  # kg/s
    m_WASSER = 60E-3
    cp_WASSER = 4.1819  # kJ/kg*K
    m_OEL_2 = 250E-3
    m_OEL_3 = 60E-3

    h_g = CP.PropsSI('H', 'P', 101325, 'Q', 1, fluid)
    h_liq = CP.PropsSI('H', 'P', 101325, 'Q', 0, fluid)
    h_v = h_g - h_liq # Vedampfungsenthalpie
    Tcrit = CP.PropsSI('TCRIT',fluid)

    """
    Pumpe: Zustand 1 so anpassen, dass Fluid bei gewünschtem Druck unterkühlt vorliegt
    """
    p1 = 100000  # Pa
    T1 = 229 # Kelvin
    etaP = 0.9  # wird hier als konstant angesehen

    h1 = CP.PropsSI('H', 'T', T1, 'P', p1, fluid)
    s1 = CP.PropsSI('S', 'T', T1, 'P', p1, fluid)
    p2 = 2000000  # Pa
    v1 = 1 / (CP.PropsSI('D', 'P', p1, 'T', T1, fluid))  # m3/kg

    w_p = (v1 * (p2 - p1)) / etaP # J
    P_p = (w_p * m_ORC)
    h2 = w_p + h1
    T2 = CP.PropsSI('T', 'H', h2, 'P', p2, fluid)
    s2 = CP.PropsSI('S', 'T', T2, 'P', p2, fluid)
    if CP.PhaseSI('H', h2, 'P', p2, fluid) == 'twophase':
        raise ValueError("no liquid state at pump outlet!")

    """ 
    Betrachtet wird ein Doppelrohrwärmeübertrager
    Innen befindet sich das Arbeitsfluid und außen das Speicherfluid
    """

    d_i = 30E-3
    d_ai = 32E-3
    d_aa = 38E-3

    i_1 = 0
    while i_1 < 3:

        '''
        Auslegung des Wärmeübertragers 1 (Unterkühlte Flüssigkeit zu siedender Flüssigkeit)
        '''
        T2_siedend = CP.PropsSI('T', 'P', p2, 'Q', 0, fluid)
        h2_siedend = CP.PropsSI('H', 'P', p2, 'Q', 0, fluid)
        s2_siedend = CP.PropsSI('S', 'P', p2, 'Q', 0, fluid)
        Q_zu1 = m_ORC * (h2_siedend - h2)
        speicherfluid1 = "REFPROP::WATER"

        Tlow_H = T2_siedend + T_pinch
        p_Tank1 = 100000  # Pa

        lambda_Wasser = 0.6
        alpha_a_1 = alpha_outside_tube(d_ai,d_aa,lambda_Wasser)
        alpha_i_1 = alpha_1P_i(p2, T2, fluid, m_ORC, d_i)

        dTB_1 = Tlow_H - T2_siedend

        lambda_Kupfer = 401 #Quelle VDI Wärmeatlas
        A_i = np.pi * d_i * l1
        A_a = np.pi * d_ai * l1
        R_konv_innen1 = 1 / (A_i * alpha_i_1)
        R_konv_aussen1 = 1 / (A_a * alpha_a_1)
        R_waermeleitung1 = np.log(d_ai / d_i) / (2 * np.pi * l1 * lambda_Kupfer)
        R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
        Tlow_L = fsolve(solver_for_WU1, T2 + 0.01, args=(Q_zu1, R_ges1, T2, dTB_1))
        m_WASSER = (Q_zu1 / 1000) / (cp_WASSER * (Tlow_H - Tlow_L)) * 1000
        i_1 += 1

    i_2 = 0
    while i_2 < 3:

        '''
        Auslegung des Wärmeübertragers 2 (siedende Flüssigkeit zu Sattdampf)
        isotherme Zustandsänderung, daher über 1.HS
        '''
        #arbeitsfluid2 = shell heat transfer oil s2
        Tmittel_L = T2_siedend + T_pinch  # K
        p_Tank2 = 100000  # Pa
        lambda_Oel = 0.129 # also konstant angenommen
        T2_sattdampf = CP.PropsSI('T', 'P', p2, 'Q', 1, fluid)
        viscosity2_liq = CP.PropsSI('VISCOSITY', 'Q', 0, 'P', p2, fluid)
        viscosity2_gas = CP.PropsSI('VISCOSITY', 'Q', 1, 'P', p2, fluid)
        cp2_liq = CP.PropsSI('C', 'P', p2, 'Q', 0, fluid)
        sigma = CP.PropsSI('SURFACE_TENSION', 'P', p2, 'Q', 0, fluid)

        h2_sattdampf = CP.PropsSI('H', 'P', p2, 'Q', 1, fluid)
        s2_sattdampf = CP.PropsSI('S', 'P', p2, 'Q', 1, fluid)
        Q_zu2 = m_ORC * (h2_sattdampf - h2_siedend)

        rho2_siedend = CP.PropsSI('D', 'Q', 0, 'P', p2, fluid)
        rho2_sattdampf = CP.PropsSI('D', 'Q', 1, 'P', p2, fluid)
        lambda_fluid_2 = CP.PropsSI('CONDUCTIVITY', 'Q', 0, 'P', p2, fluid)
        Te = T2_siedend
        dPsat = 0
        alpha_i_zweiphasig = alpha_boiling(m_ORC, 0.5, d_i, rho2_siedend, rho2_sattdampf, viscosity2_liq, viscosity2_gas, lambda_fluid_2, cp2_liq, h_v, sigma, dPsat, Te)
        alpha_a_2 = alpha_outside_tube(d_ai, d_aa, lambda_Oel)

        cp_Oel = 2.2 # in kJ/kg, konstant angenommen

        A_i_2 = np.pi * d_i * l2
        A_a_2 = np.pi * d_ai * l2

        dTA_2 = Tmittel_L - T2_siedend

        R_konv_innen2 = 1 / (A_i_2 * alpha_i_zweiphasig)
        R_konv_aussen2 = 1 / (A_a_2 * alpha_a_2)
        R_waermeleitung2 = np.log(d_ai / d_i) / (2 * np.pi * l2 * lambda_Kupfer)
        R_ges2 = R_konv_innen2 + R_konv_aussen2 + R_waermeleitung2
        Tmittel_H = fsolve(solver_for_WU2, T2_sattdampf + 0.01, args=(Q_zu2, R_ges2, T2_sattdampf, dTA_2))
        m_OEL_2 = (Q_zu2 / 1000) / (cp_Oel * (Tmittel_H - Tmittel_L)) * 1000
        i_2 += 1

    i_3 = 0
    while i_3 < 3:

        '''
        Auslegung des Wärmeübertragers 3 (Sattdampf zu überhitzten Dampf)
        '''
        T3 = 400 # K
        Thoch_H = T3 + T_pinch  # K

        dTB_3 = Thoch_H - T3
        if dTB_3 < 0:
            raise ValueError("temperature of working fluid is higher then storage fluid")

        h3 = CP.PropsSI('H', 'T', T3, 'P', p2, fluid)
        Q_zu3 = m_ORC * (h3 - h2_sattdampf) #W

        alpha_i_3 = alpha_1P_i(p2, T2_sattdampf, fluid, m_ORC, d_i)
        alpha_a_3 = alpha_outside_tube(d_ai, d_aa, lambda_Oel)

        A_i_3 = np.pi * d_i * l3
        A_a_3 = np.pi * d_ai * l3
        R_konv_innen3 = 1 / (A_i_3 * alpha_i_3)
        R_konv_aussen3 = 1 / (A_a_3 * alpha_a_3)
        R_waermeleitung3 = np.log(d_ai / d_i) / (2 * np.pi * l3 * lambda_Kupfer)
        R_ges3 = R_konv_innen3 + R_konv_aussen3 + R_waermeleitung3

        Thoch_L = fsolve(solver_for_WU3, T2_sattdampf + 0.01, args=(Q_zu3, R_ges3, T2_sattdampf, dTB_3))
        m_OEL_3 = (Q_zu3 / 1000) / (cp_Oel * (Thoch_H - Thoch_L)) * 1000
        i_3 += 1

    '''
    Überprüfung Strömungsgeschwindigkeit
    '''
    rho3 = CP.PropsSI('D', 'T', T3, 'P', p2, fluid)
    v_3 = m_ORC / rho3
    c_3 = v_3 / A_i_3
    if c_3 > 2:
        raise ValueError("the flow rate of the fluid is too high")

    Q_zu_ges = Q_zu1 + Q_zu2 + Q_zu3

    '''
    Berechnung Expander
    '''
    eta_Expander = isentroper_Wirkungsgrad()
    s3 = CP.PropsSI('S', 'P', p2, 'H', h3, fluid)
    p4 = p1
    h4s = CP.PropsSI('H', 'S', s3, 'P', p4, fluid)
    T4s = CP.PropsSI('T', 'P', p4, 'H', h4s, fluid)
    w_ts = (h4s - h3)
    w_t = w_ts * eta_Expander
    P_t = m_ORC * w_t
    h4 = eta_Expander * w_ts + h3
    T4 = CP.PropsSI("T", "H", h4, "P", p1, fluid)
    s4 = CP.PropsSI("S", "H", h4, "P", p1, fluid)

    h4_siedend = CP.PropsSI('H', 'P', p4, 'Q', 0, fluid)
    T4_siedend = CP.PropsSI('T', 'P', p4, 'H', h4_siedend, fluid)
    s4_siedend = CP.PropsSI('S', 'P', p4, 'H', h4_siedend, fluid)
    h4_sattdampf = CP.PropsSI('H', 'P', p4, 'Q', 1, fluid)
    T4_sattdampf = CP.PropsSI('T', 'P', p4, 'H', h4_sattdampf, fluid)
    s4_sattdampf = CP.PropsSI('S', 'P', p4, 'H', h4_sattdampf, fluid)

    i_k = 0
    while i_k < 1:

        '''
        Kondensator 1, ÜD -> SD, Kühlmedium R23
        '''
        kuehlmittel1 = "REFPROP::R23"
        p_Kuehlmittel1 = 1500000  # Pa

        Q_ab1 = m_ORC * (h4 - h1)
        Ta_kuehlmittel1 = T4 - 75
        ha_kuehlmittel1 = CP.PropsSI('H','T',Ta_kuehlmittel1,'P',p_Kuehlmittel1,kuehlmittel1)
        dTA_k1 = T4 - Ta_kuehlmittel1

        alpha_i_k1 = alpha_1P_i(p4,T4,fluid,m_ORC,d_i)
        alpha_a_k1 = alpha_1P_annulus(p4,Ta_kuehlmittel1,kuehlmittel1,m_Kuehlmittel1,d_ai,d_aa)

        A_i_k1 = np.pi * d_i * l_k1
        A_a_k1 = np.pi * d_ai * l_k1
        R_konv_innen_k1 = 1 / (A_i_k1 * alpha_i_k1)
        R_konv_aussen_k1 = 1 / (A_a_k1 * alpha_a_k1)
        R_waermeleitung_k1 = np.log(d_ai / d_i) / (2 * np.pi * l_k1 * lambda_Kupfer)
        R_ges_k1 = R_konv_innen_k1 + R_konv_aussen_k1 + R_waermeleitung_k1

        Te_kuehlmittel1 = fsolve(solveT_K, T1 - 1, args=(Q_ab1, R_ges_k1, T1, dTA_k1))
        he_kuehlmittel1 = CP.PropsSI('H', 'P', p_Kuehlmittel1, 'T', Te_kuehlmittel1, kuehlmittel1)
        m_Kuehlmittel1 = (Q_ab1 / (ha_kuehlmittel1 - he_kuehlmittel1)) * 1000
        i_k += 1

    Q_ab_ges = Q_ab1

    "Berechnung thermischer Wirkungsgrad"
    P_netto = abs(P_t + P_p)
    eta_th = P_netto / Q_zu_ges

    "Entropieberechnung"
    Tm1 = (Tlow_H + Tlow_L)/2
    Tm2 = (Tmittel_H + Tmittel_L)/2
    Tm3 = (Thoch_H + Thoch_L)/2
    Tmk1 = (Ta_kuehlmittel1 + Te_kuehlmittel1)/2
    s_evap1 = -Q_zu1 / Tm1
    s_evap2 = -Q_zu2 / Tm2
    s_evap3 = -Q_zu3 / Tm3
    s_kond = Q_ab_ges / Tmk1
    s_irr = s_evap1 + s_evap2 + s_evap3 + s_kond

    h_g1 = CP.PropsSI('H', 'P', p2, 'Q', 1, fluid)
    h_liq1 = CP.PropsSI('H', 'P', p2, 'Q', 0, fluid)
    h_v1 = h_g1 - h_liq1  # Vedampfungsenthalpie

    if s_irr < 0:
        raise ValueError('Die irreversibel erzeugte Entropie ist kleiner null!')

    a.append(eta_th)
    b.append(T2_sattdampf-T2_siedend)
    c.append(P_netto)
    d.append(round(x_a,1))
    f.append(s_irr)
    g.append(p2)
    h.append(T3)
    j.append(s_evap2[0])
    k.append(s_evap3[0])
    l.append(s_kond[0])
    m.append(s_evap1[0])
    s.append(Q_zu2)
    t.append(Tmittel_H-Tmittel_L)
    u.append(h_v1/1000)

plt.figure(7)
plt.plot(d,b,color='blue')
plt.title(f"Temperaturgleit für Molenbruch Isobutan\nfür m_ORC = {m_ORC*1000}g/s", fontsize=18)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=16)
plt.ylabel('Temperaturgleit [K]', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.show()

plt.figure(8)
plt.plot(d,s,color='blue')
plt.title(f"Q_zu2 für Molenbruch Isobutan\nfür m_ORC = {m_ORC*1000}g/s", fontsize=18)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=13)
plt.ylabel('Q_zu2 [W]', fontsize=13)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.grid(True)
plt.show()

plt.figure(9)
plt.plot(d,t,color='blue')
plt.title(f"delta_T_mittel für Molenbruch Isobutan\nfür m_ORC = {m_ORC*1000}g/s", fontsize=18)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=16)
plt.ylabel('Temperaturdifferenz [K]', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.show()

plt.figure(10)
plt.plot(d,u,color='blue')
plt.title(f"Verdampfungsenthalpie für Molenbruch Isobutan\nfür m_ORC = {m_ORC*1000}g/s", fontsize=18)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=16)
plt.ylabel('Verdampfungsenthalpie [kJ/kg]', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.show()

plt.figure(4)
temperature = (d)
entropy = {
    'S_evap1': m,
    'S_evap2': j,
    'S_evap3': k,
    'S_kond': l,
}
xlabel = np.arange(len(d))
barWidth = 0.25
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in entropy.items():
    offset = barWidth * multiplier
    rects = ax.bar(xlabel + offset, measurement, barWidth, label=attribute)
    multiplier +=1

ax.set_xlabel('Molenbruch Isobutan x_a',fontsize=18)
ax.set_ylabel('Entropieströme [W/K]', fontsize=18)
plt.yticks(fontsize=12)
ax.set_title('Gesamtentropieerzeugung (Sirr)\nüber Molenbruch Isobutan',fontsize=18)
ax.set_xticks(xlabel + barWidth, temperature,fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.legend(loc='upper left')
plt.show()

plt.figure(3)
plt.plot(d,a,color='blue')
plt.title(f"Thermischer Wirkungsgrad über Molenbruch Isobutan\nfür m_ORC = {m_ORC*1000}g/s", fontsize=16)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=14)
plt.ylabel('Thermischer Wirkungsgrad', fontsize=14)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.grid(True)
plt.show()

plt.figure(5)
plt.plot(d,f,color='blue')
plt.title(f"Gesamtentropieerzeugung (Sirr) über Molenbruch\nIsobutanfür m_ORC = {m_ORC*1000}g/s", fontsize=16)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=16)
plt.ylabel('Gesamtentropiestrom [W/K] ', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.show()

plt.figure(6)
plt.plot(d,c,color='blue')
plt.title(f"Nettoleistung über Molenbruch Isobutan\n für m_ORC = {m_ORC*1000}g/s", fontsize=17)
plt.xlabel('Molenbruch Isobutan x_a', fontsize=16)
plt.ylabel('Nettoleistung [W]', fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.show()

point_label = ["1", "2", "2a", "2b", "3", "4", "4a"]
x = [s1, s2, s2_siedend, s2_sattdampf, s3, s4, s4_siedend]
y = [T1, T2, T2_siedend, T2_sattdampf, T3, T4, T4_siedend]
plt.figure(1)
for i in range(len(x)):
    plt.plot(x[i], y[i], '*', markersize=15)
    plt.annotate(point_label[i], (x[i]+25, y[i]), fontsize=12)
plt.xlabel("s in J/kg/K")
plt.ylabel("T in K")
plt.legend()

# Berechnung des Nassdampfbereichs #
s_i = []
s_j = []
t_step = np.linspace(200, Tcrit-0.01, 50)
for t_i in t_step:
    s_i1 = CP.PropsSI('S', 'T', t_i, 'Q', 0, fluid)
    s_i2 = CP.PropsSI('S', 'T', t_i, 'Q', 1, fluid)
    s_i.append(s_i1)
    s_j.append(s_i2)

plt.plot(s_i, t_step, 'k-')
plt.plot(s_j, t_step, 'k-', label="wet steam region")
plt.title('T-s-Diagramm für ' + fluid)

# Berechnung Isobare #
s_step = np.linspace(200, 2700, 100)
for px in [p1, p2]:
    t_isobar = []
    for si in s_step:
        t_iso = CP.PropsSI('T', 'S', si, 'P', px, fluid)
        t_isobar.append(t_iso)

    plt.plot(s_step, t_isobar, 'b:', label="isobare")
plt.legend()
plt.show()

# h_dot-T diagram
point_label = ["1", "2", "2a", "2b", "3", "4"]
y = [T1, T2, T2_siedend, T2_sattdampf, T3, T4]
x2 = np.array([h1, h2, h2_siedend, h2_sattdampf, h3, h4])
x2 = m_ORC * x2
plt.figure(2)
for i in range(len(x2)):
    plt.subplot(111)
    plt.plot(x2[i], y[i], '*', markersize=15)
    plt.annotate(point_label[i], (x2[i]+25, y[i]), fontsize=12)
    plt.xlabel("H in W",fontsize=16)
    plt.ylabel("T in K",fontsize=16)
    plt.legend()

# Berechnung des Nassdampfbereichs #
h_i = []
h_j = []
t_step = np.linspace(200, Tcrit-0.01, 50)
for t_i in t_step:
    h_i1 = CP.PropsSI('H', 'T', t_i, 'Q', 0, fluid)
    h_i2 = CP.PropsSI('H', 'T', t_i, 'Q', 1, fluid)
    h_i.append(h_i1)
    h_j.append(h_i2)

plt.plot(np.array(h_i) * m_ORC, t_step, 'k-')
plt.plot(np.array(h_j) * m_ORC, t_step, 'k-', label="Nassdampfgebiet")
plt.title('T-H-Diagramm für\n' + fluid,fontsize=16)

# Berechnung Isobare #
h_step = np.linspace(0, 900000, 100)
for px in [p1, p2]:
    t_isobar = []
    for hi in h_step:
        t_iso = CP.PropsSI('T', 'H', hi, 'P', px, fluid)
        t_isobar.append(t_iso)

    plt.plot(h_step * m_ORC, t_isobar, '0.7', linestyle=":")
plt.legend()

# adding secondary fluids to plot figure 2

x_sec_evap = np.linspace(h1 * m_ORC, h4 * m_ORC, 100)
y_sec_evap = np.linspace(Te_kuehlmittel1, Ta_kuehlmittel1, 100)
plt.plot(x_sec_evap, y_sec_evap, 'm', label="Kondensation")
plt.legend()

x_sec_sc = np.linspace(h3 * m_ORC, (h3 + (h2_sattdampf - h3)) * m_ORC, 100)
y_sec_sc = np.linspace(Thoch_H, Thoch_L, 100)
plt.plot(x_sec_sc, y_sec_sc, 'r')
plt.plot(x_sec_sc, y_sec_sc, 'r',label="Wärmeübertrager 3")

x_sec_ws = np.linspace(h2_sattdampf * m_ORC, (h2_sattdampf + (h2_siedend - h2_sattdampf)) * m_ORC, 100)
y_sec_ws = np.linspace(Tmittel_H, Tmittel_L, 100)
plt.plot(x_sec_ws, y_sec_ws, 'g')
plt.plot(x_sec_ws, y_sec_ws, 'g',label="Wärmeübertrager 2")

x_sec_sh = np.linspace(h2_siedend * m_ORC, (h2_siedend + (h2 - h2_siedend)) * m_ORC, 100)
y_sec_sh = np.linspace(Tlow_H, Tlow_L, 100)
plt.plot(x_sec_sh, y_sec_sh, 'b')
plt.plot(x_sec_sh, y_sec_sh, 'b',label="Wärmeübertrager 1")

plt.legend()
ax.legend(loc='upper left')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

T_L = (Ta_kuehlmittel1 - Te_kuehlmittel1) / np.log(Ta_kuehlmittel1/Te_kuehlmittel1)
T_H = (Thoch_H - Thoch_L) / np.log(Thoch_H/Thoch_L)
eta_C = 1 - T_L/T_H