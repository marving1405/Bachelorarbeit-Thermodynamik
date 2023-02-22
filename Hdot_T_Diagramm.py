# h_dot-T diagram
point_label = ["1", "2", "2siedend", "2sattdampf", "3", "4", "4siedend"]
y = [T1, T2, T2_siedend, T2_sattdampf, T3, T4, T4_siedend]
x2 = np.array([h1, h2, h2_siedend, h2_sattdampf, h3, h4, h4_siedend])
x2 = m_ORC * x2
plt.figure(2)
for i in range(len(x2)):
    plt.plot(x2[i], y[i], '*', markersize=15)
    plt.annotate(point_label[i], (x2[i]+25, y[i]), fontsize=12)
plt.xlabel("h_dot in J/s")
plt.ylabel("T in K")
plt.legend()

# Berechnung des Nassdampfbereichs #
h_i = []
h_j = []
t_step = np.linspace(250, 442, 50)
for t_i in t_step:
    h_i1 = CP.PropsSI('H', 'T', t_i, 'Q', 0, fluid)
    h_i2 = CP.PropsSI('H', 'T', t_i, 'Q', 1, fluid)
    h_i.append(h_i1)
    h_j.append(h_i2)

plt.plot(np.array(h_i) * m_ORC, t_step, 'k-')
plt.plot(np.array(h_j) * m_ORC, t_step, 'k-', label="wet steam region")
#plt.xlabel('h in J/kg')
#plt.ylabel('T in K')
plt.title('T-h-Diagramm für ' + fluid)

# Berechnung Isobare #
h_step = np.linspace(0, 700000, 100)
for px in [p1, p2]:
    t_isobar = []
    for hi in h_step:
        t_iso = CP.PropsSI('T', 'H', hi, 'P', px, fluid)
        t_isobar.append(t_iso)

    plt.plot(h_step * m_dot, t_isobar, 'b:', label="isobare")
plt.legend()
