# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   OPT FUNCTION INITIALIZATION         -------------------------------------> #
fluids = ["Water", "Carbon Dioxide", "Ethane", "Methane", "Argon", "Neon", "Helium"]
acntr_factors = list()

for fluid in fluids:

    tp_vap = PlantThermoPoint([fluid], [1], unit_system="MASS BASE SI")
    tp_vap.set_variable("P", tp_vap.RPHandler.PC)
    tp_vap.set_variable("T", tp_vap.RPHandler.TC)
    acntr_factors.append(tp_vap.evaluate_RP_code("ACF"))

p_list = np.logspace(-3, 0, 40000)


# %%-------------------------------------   EVALUATE RP VALUES                  -------------------------------------> #
pbar = tqdm(desc="Calculating Points", total=len(fluids)*len(p_list))
t_sat_arr = np.empty(np.array([len(fluids), len(p_list)]))
t_sat_arr[:, :] = np.nan

for i in range(len(fluids)):

    tp_vap = PlantThermoPoint([fluids[i]], [1], unit_system="MASS BASE SI")
    tp_liq = tp_vap.duplicate()

    p_crit = tp_vap.RPHandler.PC
    t_crit = tp_vap.RPHandler.TC

    for j in range(len(p_list)):

        tp_vap.set_variable("P", p_list[j]*p_crit)
        tp_vap.set_variable("Q", 1)
        tp_liq.set_variable("P", p_list[j] * p_crit)
        tp_liq.set_variable("Q", 0)

        if tp_vap.get_variable("T") / t_crit > 0.4:
            t_sat_arr[i, j] = tp_vap.get_variable("T") / t_crit

        pbar.update(1)

pbar.close()


# %%-------------------------------------   OPTIMIZE RATIO                      -------------------------------------> #
def opt_function(i_curr, minimise_max=False):

    def funct(x):

        a = x[0]
        b = np.log(0.7) * a / (np.power(10, - a * (acntr_factors[i_curr] + 1)) - 1)
        t_sat_actr = np.exp(b / a * (np.power(p_list, a) - 1))

        if minimise_max:
            return np.nanmax(np.power(t_sat_actr - t_sat_arr[i_curr, :], 2) / t_sat_arr[i_curr, :])

        return np.nanmean(np.power(t_sat_actr - t_sat_arr[i_curr, :], 2) / t_sat_arr[i_curr, :])

    return funct

ratios = list()
coeff_list = list()
t_sat_approx = np.empty(np.array([len(fluids), len(p_list)]))
t_sat_approx[:, :] = np.nan

pbar = tqdm(desc="Optimizing Ratios", total=len(fluids))
for i in range(len(fluids)):

    res = minimize(opt_function(i), np.array([0.11]))
    ratios.append(res.x[0])
    a = res.x[0]
    b = np.log(0.7) * a / (np.power(10, - a * (acntr_factors[i] + 1)) - 1)

    t_sat_approx[i, :] = np.exp(b / a * (np.power(p_list, a) - 1))
    coeff_list.append({"a": a, "b": b})

    pbar.update(1)

pbar.close()
rel_error = (t_sat_approx - t_sat_arr) / t_sat_arr


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
ax_err = ax.twinx()

for i in range(len(fluids)-1):

    line = ax.plot(t_sat_arr[i, :], p_list, label=fluids[i])[0]
    ax.plot(t_sat_approx[i, :], p_list, "--", color=line.get_color(), alpha=0.75)
    ax_err.plot(t_sat_approx[i, :], rel_error[i, :] * 100, "-.", color=line.get_color(), alpha=0.50)

ax.set_title("Saturation Line - Approximation")
ax.set_xlabel("$T_{rel}$ [-]")
ax.set_ylabel("$p_{rel}$ [-]")
ax_err.set_ylabel("Error $T_{sat}}$ [%]")
ax.set_yscale("log")
ax.legend()
plt.show()

# %%-------------------------------------   PLOT RATIOS                         -------------------------------------> #
plt.scatter(acntr_factors, ratios)
plt.show()
