Kondensator
2
Nassdampfbereich

l_k2 = 45  # m
kuehlmittel2 = 'REFPROP:METHANOL'
p_kuehlmittel2 = p_Kuehlmittel1
m_Kuehlmittel2 = 100E-3
h4_siedend = CP.PropsSI('H', 'P', p4, 'Q', 0, fluid)
T4_siedend = CP.PropsSI('T', 'P', p4, 'H', h4_siedend, fluid)
s4_siedend = CP.PropsSI('S', 'P', p4, 'H', h4_siedend, fluid)
Q_ab2 = m_ORC * (h4_sattdampf - h1)
Ta_kuehlmittel2 = T4_sattdampf - 5  # K
ha_kuehlmittel2 = CP.PropsSI('H', 'T', Ta_kuehlmittel2, 'P', p4, kuehlmittel2)

dTA_k2 = T4_sattdampf - Ta_kuehlmittel2

rhol = CP.PropsSI('D', 'Q', 0, 'P', p4, fluid)
rhog = CP.PropsSI('D', 'Q', 1, 'P', p4, fluid)
lambda_liq = CP.PropsSI('CONDUCTIVITY', 'Q', 0, 'P', p_kuehlmittel2, kuehlmittel2)
viscosity_liq = CP.PropsSI('VISCOSITY', 'Q', 0, 'P', p4, fluid)
cp_liq = CP.PropsSI('C', 'Q', 0, 'P', p4, kuehlmittel2)

alpha_i_nassdampf = alpha_condensation(m_ORC, rhog, rhol, lambda_liq, viscosity_liq, cp_liq, d_i, 0.5)
alpha_a_k2 = alpha_1P_annulus(p4, Ta_kuehlmittel2, kuehlmittel2, m_Kuehlmittel2, d_ai, d_aa)

A_i_k2 = np.pi * d_i * l_k2
A_a_k2 = np.pi * d_ai * l_k2
R_konv_innen_k2 = 1 / (A_i_k2 * alpha_i_nassdampf)
R_konv_aussen_k2 = 1 / (A_a_k2 * alpha_a_k2)
R_waermeleitung_k2 = np.log(d_ai / d_i) / (2 * np.pi * l_k2 * lambda_Kupfer)
R_ges_k2 = R_konv_innen_k2 + R_konv_aussen_k2 + R_waermeleitung_k2

Te_kuehlmittel2 = fsolve(solveT_K, Ta_kuehlmittel2 - 10, args=(Q_ab2, R_ges_k2, T1, dTA_k2))
he_kuehlmittel2 = CP.PropsSI('H', 'P', p_kuehlmittel2, 'T', Te_kuehlmittel2, kuehlmittel2)