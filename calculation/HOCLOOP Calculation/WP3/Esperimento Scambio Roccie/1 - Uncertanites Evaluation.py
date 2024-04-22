# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np


# %%------------   INIT RANDOM ARRAYS                     -----------------------------------------------------------> #
n_points = 10000000

alpha_mu = -6
alpha_sd = 1
exp_alpha = np.random.normal(alpha_mu, alpha_sd, n_points)
alpha_rand = np.power(10, exp_alpha)

k_mu = -10
k_sd = 2
exp_k = np.random.normal(k_mu, k_sd, n_points)
k_rand = np.power(10, exp_k)


# %%------------   EVALUATE Pe_MAX                        -----------------------------------------------------------> #
r_list = np.array([0.001, 0.01, 0.1])
Pe_max_list_ovr = list()
t_d_pe_list_ovr = list()
t_d_perc_list_ovr = list()

dpdl = 10 ** 4
mu = 10 ** -3
a = 12.6357261129
b = -1.5360002093
Pe_perc_list = [0.25, 0.5, 0.75, 1]

for Pe_perc in Pe_perc_list:

    print(Pe_perc)

    Pe_max_list = list()
    t_d_pe_list = list()
    t_d_perc_list = list()

    for r in r_list:

        Pe_max = 2 * r * k_rand * dpdl / (mu * alpha_rand)
        t_d_min = alpha_rand / (r ** 2)
        t_d_max = t_d_min * 10 ** 4
        t_d_pe_max = a * np.power(Pe_max*Pe_perc, b)
        t_d_pe_max_perc = np.log10(t_d_pe_max / t_d_min) / np.log10(t_d_max / t_d_min)
        print(np.count_nonzero((t_d_pe_max_perc < 1) & (t_d_pe_max_perc > 0)) / len(t_d_max))

        Pe_max_list.append(Pe_max)
        t_d_pe_list.append([t_d_min, t_d_pe_max, t_d_max])
        t_d_perc_list.append(t_d_pe_max_perc)

    Pe_max_list_ovr.append(Pe_max_list)
    t_d_pe_list_ovr.append(t_d_pe_list)
    t_d_perc_list_ovr.append(t_d_perc_list)


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
plot_linear = False
for i in range(len(r_list)):

    plt_array = t_d_pe_list_ovr[-1][0][i]
    hist, logbins = np.histogram(plt_array, bins=300, density=True)

    if not plot_linear:
        logbins = np.logspace(np.log10(logbins[0]),np.log10(logbins[-1]),len(logbins))
        plt.xscale("log")

    plt.hist(plt_array, bins=logbins, alpha=0.5, label="$r_0$={} cm".format(r_list[i]*100))

plt.legend()
plt.show()


# %%------------   PLOT T_D_%_list                        -----------------------------------------------------------> #
plot_linear = True
plt.figure(figsize=(10,4))
plt.subplot(121)
for i in range(len(r_list)):

    plt_array = t_d_perc_list_ovr[3][i]
    if plot_linear:
        hist, logbins = np.histogram(plt_array, bins=np.linspace(0, 1, 500), density=True)
    else:
        plt.xscale("log")
        hist, logbins = np.histogram(plt_array, bins=np.logspace(-3, 0, 150), density=True)

    yhat = savgol_filter(hist, 100, 2)
    plt.plot(logbins[:-1], yhat, alpha=0.5, label="$r_0$={} cm".format(r_list[i]*100))

plt.xlabel("${t_d}_{pos\%}$")
plt.ylabel("probability density [-]")
plt.title("$r_0$ effect")
plt.legend()

plt.subplot(122)
for i in range(len(Pe_perc_list)):

    plt_array = t_d_perc_list_ovr[i][0]
    if plot_linear:
        hist, logbins = np.histogram(plt_array, bins=np.linspace(0, 1, 500), density=True)
    else:
        plt.xscale("log")
        hist, logbins = np.histogram(plt_array, bins=np.logspace(-3, 0, 150), density=True)

    yhat = savgol_filter(hist, 100, 2)
    plt.plot(logbins[:-1], yhat, alpha=0.5, label="$Pe_\%$={}%".format(Pe_perc_list[i]*100))

plt.xlabel("${t_d}_{pos\%}$")
ax = plt.gca()
ax.get_yaxis().set_visible(False)
plt.title("$Pe_\%$ effect")
plt.legend()

plt.tight_layout(pad=2)
plt.show()

