# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   EVALUATE RP VALUES                  -------------------------------------> #
fluids = ["Water", "Carbon Dioxide", "Ethane", "Methane", "Argon", "Neon", "Helium"]
acntr_factors = [0.344, 0.225, 0.0995, 0.001142, -0.00219, -0.0355, -0.3836]

p_list = np.logspace(-3, 0, 40000)

t_sat_arr = np.empty(np.array([len(fluids), len(p_list)]))
dh_vap_arr = np.empty(np.array([len(fluids), len(p_list)]))
dv_vap_arr = np.empty(np.array([len(fluids), len(p_list)]))
ratio_arr = np.empty(np.array([len(fluids), len(p_list)]))

t_sat_arr[:, :] = np.nan
dh_vap_arr[:, :] = np.nan
dv_vap_arr[:, :] = np.nan
ratio_arr[:, :] = np.nan

pbar = tqdm(desc="Calculating Points", total=len(fluids)*len(p_list))
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
            dh_vap_arr[i, j] = tp_vap.get_variable("H") - tp_liq.get_variable("H")
            dv_vap_arr[i, j] = 1 / tp_vap.get_variable("rho") - 1/ tp_liq.get_variable("rho")

            if j < len(p_list) - 1:

                ratio_arr[i, j] = dv_vap_arr[i, j] / dh_vap_arr[i, j]

            else:

                tp_vap.set_variable("P", p_crit)
                tp_vap.set_variable("T", t_crit)

                rho = tp_vap.get_variable("rho")
                drho_dp = tp_vap.get_derivative("rho", "P", "T")
                drho_dt = tp_vap.get_derivative("rho", "T", "P")
                d2rho_dp2 = tp_vap.get_second_derivative("rho", "P")
                d2rho_dpdt = tp_vap.get_second_derivative("rho", "P", "T")

                dv_dp = - drho_dp / (rho ** 2)
                dv_dt = - drho_dt / (rho ** 2)
                d2v_dp2 = drho_dp ** 2 / (2 * rho ** 3) - d2rho_dp2 / (rho ** 2)
                d2v_dpdt = drho_dp * drho_dt / (2 * rho ** 3) - d2rho_dpdt / (rho ** 2)
                d2h_dp2 = dv_dp - t_crit * d2v_dpdt

                ratio_arr[i, j] = d2v_dp2 / d2h_dp2
                # ratio_arr[i, j] = - drho_dp / (drho_dt * t_crit)

        pbar.update(1)

pbar.close()


# %%-------------------------------------   CHECK LINEAR REGRESSION             -------------------------------------> #
models = list()
plot_graph = False
approx_index = int(len(p_list) - 1)

t_sat_approx = np.empty(np.array([len(fluids), len(p_list)]))
t_sat_approx[:, :] = np.nan

x_values = np.log(p_list).reshape(-1, 1)[:-1]

for i in range(len(fluids)):

    tp = PlantThermoPoint([fluids[i]], [1], unit_system="MASS BASE SI")
    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    # Interpolate Values
    y_values = np.log(ratio_arr[i, :-1] / 1)
    model = LinearRegression(fit_intercept=True)
    exclude_values = np.logical_not(np.isnan(y_values))

    model.fit(

        x_values[exclude_values],
        y_values[exclude_values]

    )

    if plot_graph:
        plt.title(fluids[i])
        plt.plot(x_values, y_values)
        plt.plot(x_values, model.predict(x_values))
        plt.show()

    models.append(model)

    # Integrate Saturation Temperatures
    a = model.coef_[0] + 1
    dvdh_crit = np.exp(models[i].intercept_)

    t_ref = t_sat_arr[i, approx_index]
    p_ref = p_list[approx_index]
    t_sat_approx[i, :] = t_ref * np.exp(dvdh_crit * p_crit / a * (np.power(p_list, a) - np.power(p_ref, a)))

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


# %%-------------------------------------   REGRESSION MODELS WITH ACENTRICITY  -------------------------------------> #
fig, ax_ratio = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
ax_inter = ax_ratio.twinx()
ratios = list()
intercepts = list()
dcp0s = list()

for i in range(len(fluids)):

    ratios.append(models[i].coef_[0])
    intercepts.append(models[i].intercept_)

    tp = PlantThermoPoint([fluids[i]], [1], unit_system="MASS BASE SI")
    tp.set_variable("P", 1e5)

    tp.set_variable("T", tp.RPHandler.TC - 0.1)
    cp_0 = tp.evaluate_RP_code("CP0")

    tp.set_variable("T", tp.RPHandler.TC + 0.1)
    dcp0s.append((tp.evaluate_RP_code("CP0") - cp_0) / 0.2)

x_values = np.concatenate((acntr_factors, dcp0s)).reshape((2, len(acntr_factors))).T
y_values = np.array(intercepts).reshape(-1, 1)
model_inter = LinearRegression(fit_intercept=True)
model_inter.fit(

    x_values,
    y_values

)

inter_errors = (model_inter.predict(x_values) - y_values) / y_values

pr = PolynomialFeatures(degree = 2)
x_values_ratio = pr.fit_transform(x_values)
y_values_ratio = np.array(ratios).reshape(-1, 1)
model_ratio = LinearRegression(fit_intercept=True)
model_ratio.fit(

    x_values_ratio,
    y_values_ratio

)
ratio_errors = (model_ratio.predict(x_values_ratio) - y_values_ratio) / y_values_ratio

ax_ratio.plot(acntr_factors, ratios, color="tab:blue")
ax_inter.plot(acntr_factors, intercepts, color="tab:orange")

acnt_test = np.concatenate((

    np.linspace(np.min(acntr_factors), np.max(acntr_factors), 50),
    np.zeros(50)

)).reshape((2, 50)).T

ax_ratio.plot(acnt_test[:,0], model_ratio.predict(pr.fit_transform(acnt_test)), "-.", color="tab:blue")
ax_inter.plot(acnt_test[:,0], model_inter.predict(acnt_test), "-.", color="tab:orange")

plt.show()


# %%-------------------------------------   TEST ACENTRICITY REGRESSION         -------------------------------------> #
plot_base = True
t_sat_actr = np.empty(np.array([len(fluids), len(p_list)]))
t_sat_actr[:, :] = np.nan
for i in range(len(fluids)):

    tp = PlantThermoPoint([fluids[i]], [1], unit_system="MASS BASE SI")
    p_crit = tp.RPHandler.PC

    # Integrate Saturation Temperatures
    if plot_base:
        a = models[i].coef_[0] + 1
        b = np.log(0.7) * a / (np.power(10, - a * (acntr_factors[i] + 1)) - 1)

    else:
        a = model_ratio.predict(pr.fit_transform([[acntr_factors[i], dcp0s[i]]])) + 1
        b = np.exp(model_inter.predict([[acntr_factors[i], dcp0s[i]]])) * p_crit

    t_sat_actr[i, :] = np.exp(b / a * (np.power(p_list, a) - 1))

rel_error_actr = (t_sat_actr - t_sat_arr) / t_sat_arr


# %%-------------------------------------   PLOT RESULTS ACENTRICITY REGRESSION -------------------------------------> #
fig, ax = plt.subplots(1, 1, dpi=300)
fig.set_size_inches(10, 5)
ax_err = ax.twinx()

for i in range(len(fluids)-1):

    line = ax.plot(t_sat_arr[i, :], p_list, label=fluids[i])[0]
    ax.plot(t_sat_actr[i, :], p_list, "--", color=line.get_color(), alpha=0.75)
    ax_err.plot(t_sat_actr[i, :], rel_error_actr[i, :] * 100, "-.", color=line.get_color(), alpha=0.50)

ax.set_title("Saturation Line - Approximation")
ax.set_xlabel("$T_{rel}$ [-]")
ax.set_ylabel("$p_{rel}$ [-]")
ax_err.set_ylabel("Error $T_{sat}}$ [%]")
ax.set_yscale("log")
ax.legend()
plt.show()