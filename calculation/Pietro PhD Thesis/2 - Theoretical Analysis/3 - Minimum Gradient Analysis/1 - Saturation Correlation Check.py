# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   EVALUATE RP VALUES                  -------------------------------------> #
fluids = ["Water", "Carbon Dioxide", "Methane"]
p_list = np.logspace(-2, 0, 500)
t_sat_arr = np.empty(np.array([len(fluids), len(p_list)]))
t_sat_arr[:, :] = np.nan

pbar = tqdm(desc="Calculating Points", total=len(fluids)*len(p_list))
for i in range(len(fluids)):

    tp = PlantThermoPoint([fluids[i]], [1], unit_system="MASS BASE SI")
    p_crit = tp.RPHandler.PC
    t_crit = tp.RPHandler.TC

    for j in range(len(p_list)):

        tp.set_variable("P", p_list[j]*p_crit)
        tp.set_variable("Q", 0)

        t_sat_arr[i, j] = tp.get_variable("T") / t_crit

        pbar.update(1)

pbar.close()


# %%-------------------------------------   EVALUATE TAYLOR EXPANSION           -------------------------------------> #
t_tay_arr = np.empty(np.array([2, len(fluids), len(p_list)]))
t_tay_arr[:, :, :] = np.nan
approx_index = int(len(p_list) / 2)
for i in range(len(fluids)):

    tp_vap = PlantThermoPoint([fluids[i]], [1], unit_system="MASS BASE SI")
    tp_liq = tp_vap.duplicate()

    p_crit = tp_vap.RPHandler.PC
    t_crit = tp_vap.RPHandler.TC

    tp_vap.set_variable("P", p_list[approx_index] * p_crit)
    tp_vap.set_variable("Q", 1)
    tp_liq.set_variable("P", p_list[approx_index] * p_crit)
    tp_liq.set_variable("Q", 0)

    p_0 = p_list[approx_index] * p_crit
    t_0 = tp_liq.get_variable("T")

    dh_vap = tp_vap.get_variable("h") - tp_liq.get_variable("h")
    dv_vap = 1 / tp_vap.get_variable("rho") - 1 / tp_liq.get_variable("rho")

    dt_dp = t_0 * dv_vap / dh_vap

    dp_rel = (p_list[:] * p_crit - p_0) / p_0
    log_dp = np.log(1 + dp_rel)

    t_tay_arr[0, i, :] = (p_0 * dt_dp * log_dp + t_0) / t_crit

    tp_vap.set_variable("P", p_list[approx_index + 1] * p_crit)
    tp_vap.set_variable("Q", 1)
    tp_liq.set_variable("P", p_list[approx_index + 1] * p_crit)
    tp_liq.set_variable("Q", 0)

    p_1 = p_list[approx_index + 1] * p_crit
    dp_curr = p_1 - p_0

    dh2_dp = (tp_vap.get_variable("h") - tp_liq.get_variable("h") - dh_vap) / dp_curr
    dv2_dp = (1 / tp_vap.get_variable("rho") - 1 / tp_liq.get_variable("rho") - dv_vap) / dp_curr
    k = (dv_vap / dh_vap + dv2_dp / dv_vap - dh2_dp / dh_vap)

    t_tay_arr[1, i, :] = (p_0 * log_dp * dt_dp * (1 + (1 + p_0 * k) * log_dp / 2) + t_0) / t_crit



# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)

for i in range(len(fluids)):

    line = ax.plot(t_sat_arr[i, :], p_list, label=fluids[i])[0]
    ax.plot(t_tay_arr[0, i, :], p_list, "--", color=line.get_color(), alpha=0.3)
    ax.plot(t_tay_arr[1, i, :], p_list, "-.", color=line.get_color(), alpha=0.6)

ax.set_title("Saturation Line - Taylor Approximation")
ax.set_xlabel("$T_{rel}$ [-]")
ax.set_ylabel("$p_{rel}$ [-]")
ax.set_yscale("log")
ax.legend()
plt.show()