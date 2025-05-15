# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   OPT FUNCTION INITIALIZATION         -------------------------------------> #
# fluids = ["Carbon Dioxide", "Water", "Propane", "R1234yf", "Ethane"]
fluids = ["Carbon Dioxide", "Ethane"]
acntr_factors = list()
gamma_in_list = list()
p_in_list = list()
c_p_list = list()
t_in = 20 + 273.15
dz = np.logspace(2, 4, 10)
dz_nd = np.zeros((len(fluids), len(dz)))
beta_preds = np.zeros((len(fluids), len(dz)))

for i in range(len(fluids)):

    fluid = fluids[i]
    tp_vap = PlantThermoPoint([fluid], [1], unit_system="MASS BASE SI")
    tp_vap.set_variable("P", tp_vap.RPHandler.PC)
    tp_vap.set_variable("T", tp_vap.RPHandler.TC)
    acntr_factors.append(tp_vap.evaluate_RP_code("ACF"))

    tp_vap.set_variable("T", t_in)
    tp_vap.set_variable("Q", 0)
    p_in = tp_vap.get_variable("P") * 1.0001

    tp_vap.set_variable("T", t_in)
    tp_vap.set_variable("P", p_in)

    rho_in = tp_vap.get_variable("rho")
    cp_in = tp_vap.get_variable("cp")
    gamma_in = rho_in * cp_in * t_in / p_in

    dz_nd[i, :] = 9.81 * dz / (cp_in * t_in)
    beta_preds[i, :] = (gamma_in * np.log(dz_nd[i, :] + 1) / (dz_nd[i, :] + 1)) + 1
    c_p_list.append(cp_in)
    gamma_in_list.append(gamma_in)
    p_in_list.append(p_in / tp_vap.RPHandler.PC)

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

phi_preds = np.zeros((len(fluids), len(dz)))
grad_t_nd_lim = np.zeros((len(fluids), len(dz)))
grad_t_lim = np.zeros((len(fluids), len(dz)))

pbar = tqdm(desc="Optimizing Ratios", total=len(fluids))
for i in range(len(fluids)):

    res = minimize(opt_function(i), np.array([0.11]))
    ratios.append(res.x[0])
    a = res.x[0]
    b = np.log(0.7) * a / (np.power(10, - a * (acntr_factors[i] + 1)) - 1)

    phi_preds[i, :] = b * p_in_list[i] ** a * (np.power(beta_preds[i, :], a) - 1)
    grad_t_nd_lim[i, :] = (np.exp(phi_preds[i, :]) - 1) / dz_nd[i, :]
    grad_t_lim[i, :] = grad_t_nd_lim[i, :] * 9.81 / c_p_list[i]

    t_sat_approx[i, :] = np.exp(b / a * (np.power(p_list, a) - 1))
    coeff_list.append({"a": a, "b": b})

    pbar.update(1)

pbar.close()
rel_error = (t_sat_approx - t_sat_arr) / t_sat_arr


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
i = 0
logy = False

ax.plot(dz, grad_t_lim[i, :]*10000)

ax.set_title("Required Gradient for {}".format(fluids[i]))
ax.set_xlabel("$Depth$ [m]")
ax.set_ylabel("$\\nabla T_{rocks_{lim}}$ [°C/km]")
ax.set_xscale("log")

if logy:
    ax.set_yscale("log")
plt.show()


# %%-------------------------------------   PLOT RATIOS                         -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
logy = False

for i in range(len(fluids)):

    ax.plot(dz, grad_t_lim[i, :]*10000, label=fluids[i])

ax.set_title("Required Gradient for different fluids")
ax.set_xlabel("$Depth$ [m]")
ax.set_ylabel("$\\nabla T_{rocks_{lim}}$ [°C/km]")
ax.set_xscale("log")
ax.legend()

if logy:
    ax.set_yscale("log")

plt.show()
