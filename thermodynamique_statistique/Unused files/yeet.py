import numpy as np
# import graphinglib as gl
# from vpython import *


# complete_p = np.loadtxt("thermodynamique_statistique/data/drude_p.csv", delimiter=",")
# mean_p = complete_p.mean(axis=1)
# t = np.arange(complete_p.shape[0])

# fig = gl.Figure(size=(11,7), x_label="Temps [u. arb.]", y_label="Quantité de mouvement [u. arb.]", figure_style="dim")
# mean_p_curve = gl.Curve(t, mean_p, line_width=2, label="Moyenne des électrons", color="yellow")
# single_p_curve = gl.Curve(t, complete_p[:,11], line_width=2, label="Électron unique", color="blue")
# fit = gl.FitFromFunction(
#     lambda t, tau, k: (mean_p[0] - k) * np.exp(-t / tau) + k,
#     mean_p_curve,
#     guesses=[10000,3.2e-25],
#     color="red"
# )
# fit.label = rf"Amortissement ajusté, $\tau={fit.parameters[0]:.2f}$"

# fig.add_elements(mean_p_curve, single_p_curve, fit)
# fig.show()

# E = vector(0,5e-23,0) # electric field
# dt = 1E-5  # pas d'incrémentation temporel
# def apply_electric_field(hitlist: list):
#     global p
#     for i in range(len(p)):
#         if i not in hitlist:
#             p[i] += E * dt


# data = {
#     magnitude : np.loadtxt(f"thermodynamique_statistique/data/drude_apos_E={magnitude}.csv", delimiter=",")
#     for magnitude in [1e-23, 3e-23, 5e-23, 7e-23, 9e-23]
# }
# t = np.arange(data[1e-23].shape[0])

# fig = gl.Figure(size=(11,7), x_label="Temps [u. arb.]", y_label="Position moyenne [u. arb.]", figure_style="dim")
# curves = []
# for (magnitude, data_array), color in zip(data.items(), ["purple", "pink", "lime", "cyan", "orange"]):
#     curves.append(gl.Curve(t, data_array[:,::2].mean(axis=1), color=color, line_style="-", label=rf"$E={magnitude}$"))
#     curves.append(gl.Curve(t, data_array[:,1::2].mean(axis=1), color=color, line_style=":"))

# fig.add_elements(*curves)
# fig.show()


targeted_particle = np.array([[ 1.00000000e-05, 7.10909945e-03, 8.68681213e-03, 0.00000000e+00],
    [ 2.00000000e-05, 1.42181989e-02, 1.73736243e-02, 0.00000000e+00],
    [ 1.00000000e-05, 7.10909945e-03, 8.68681213e-03, 0.00000000e+00],
    [ 1.00000000e-05, 7.10909945e-03, 8.68681213e-03, 0.00000000e+00],
    [ 1.00000000e-05, 7.10909945e-03, 8.68681213e-03, 0.00000000e+00],
    [ 1.00000000e-05, 7.10909945e-03, 8.68681213e-03, 0.00000000e+00],
    [ 4.00000000e-05, 3.21022448e-02, 2.60675442e-02, 0.00000000e+00],
    [ 7.00000000e-05, -7.68843323e-02, 3.72285232e-02, 0.00000000e+00],
    [ 1.00000000e-05, -1.09834760e-02, 5.31836045e-03, 0.00000000e+00],
    [ 1.00000000e-05, -1.18477422e-03, 2.68219755e-03, 0.00000000e+00],
    [ 3.00000000e-05, -3.55432267e-03, 8.04659265e-03, 0.00000000e+00]])


tau = targeted_particle[:, 0]
lmoy = np.linalg.norm(targeted_particle[:, 1:], axis=1)
print(lmoy)