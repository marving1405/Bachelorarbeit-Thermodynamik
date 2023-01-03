from calculate_alpha_aw import alpha_1P_i

def abfrage():
    abfrage = input("Rekuperator verwenden? True or False?")
    if abfrage == "True":
        alpha_R = alpha_1P_i(p4, T4, fluid, m_ORC, d_i)
    if abfrage == "False":
        return False


