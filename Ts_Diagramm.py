
point_label = ["1", "2", "2b", "2c", "3", "4", "4b", "4c"]
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
t_step = np.linspace(200, 442, 50)
for t_i in t_step:
    s_i1 = CP.PropsSI('S', 'T', t_i, 'Q', 0, fluid)
    s_i2 = CP.PropsSI('S', 'T', t_i, 'Q', 1, fluid)
    s_i.append(s_i1)
    s_j.append(s_i2)

plt.plot(s_i, t_step, 'k-')
plt.plot(s_j, t_step, 'k-', label="wet steam region")
#plt.xlabel('s in J/kg/K')
#plt.ylabel('T in K')
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