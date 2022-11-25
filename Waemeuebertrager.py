
def waermeuebertrager(l,d_i, d_ai, d_aa, lambda_fluid,Te, Ta)


    A_i = 2 * np.pi * d_i / 2 * l / 3
    A_a = 2 * np.pi * d_ai / 2 * l / 3
    R_konv_innen1 = 1 / A_i * alpha_i
    R_konv_aussen1 = 1 / A_a * alpha_a
    R_waermeleitung1 = np.log(d_aa / d_ai) / (2 * np.pi * l * lambda_fluid)
    T_siedend = PropsSI('T', 'P', p2, 'Q', 0, fluid)
    delta_T1 = T_siedend - T2
    R_ges1 = R_konv_innen1 + R_konv_aussen1 + R_waermeleitung1
    Q_zu1 = (1 / R_ges1) * delta_T1

    return Q_zu