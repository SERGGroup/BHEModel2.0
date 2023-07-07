# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.support_functions import get_np_array
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   INITIALIZATION                      -------------------------------------> #
input_point = PlantThermoPoint(["CarbonDioxide"], [1])
bhe = SimplifiedBHE(

    input_thermo_point=input_point,
    dz_well=0.1, t_rocks=30

)

n_points = 200
p_points = get_np_array(4, 50, n_points, log_scale=True)  # in MPax
h_points = get_np_array(200, 450, n_points, log_scale=False)  # in kJ/kg
res_points = get_np_array(0, 0, n_points)

h_mesh, p_mesh = np.meshgrid(h_points, p_points, indexing='ij')


# %%-------------------------------------   CALCULATION                         -------------------------------------> #
pbar = tqdm(desc="C0 map calculation", total=n_points * n_points)
res, b = np.meshgrid(res_points, h_points, indexing='ij')

for i in range(n_points):

    for j in range(n_points):

        input_point.set_variable("h", h_mesh[i, j])
        input_point.set_variable("P", p_mesh[i, j])

        bhe.update_simplified()

        res[i, j] = np.log10(bhe.C0[0])
        pbar.update(1)

pbar.close()


# %%-------------------------------------   GET SATURATION LINE                 -------------------------------------> #
support_point = PlantThermoPoint(["CarbonDioxide"], [1])

n_sat_points = 500
p_values = get_np_array(4, input_point.RPHandler.PC, n_sat_points, log_scale=True)
sat_h = get_np_array(0, 0, 2 * n_sat_points)
sat_p = get_np_array(0, 0, 2 * n_sat_points)

pbar_sat = tqdm(desc="loop", total=n_sat_points)

for i in range(n_sat_points):
    sat_p[i] = p_values[i]
    sat_p[-1 - i] = p_values[i]

    support_point.set_variable("P", p_values[i])

    support_point.set_variable("Q", 0)
    sat_h[i] = support_point.get_variable("h")

    support_point.set_variable("Q", 1)
    sat_h[-1 - i] = support_point.get_variable("h")

    pbar_sat.update(1)

pbar_sat.close()

# %%-------------------------------------   GET ISO-T LINES                     -------------------------------------> #
support_point = PlantThermoPoint(["CarbonDioxide"], [1])

n_isolines_points = 500
h_range = get_np_array(200, 450, n_isolines_points)  # in kJ/kg
T_range = get_np_array(5, 125, 25)  # in °C
P_max = 50

isolines = {

    "T": dict(),
    "rho": dict(),

}

pbar_iso_t = tqdm(desc="iso-T lines calculation", total=n_isolines_points * len(T_range))

for T in T_range:

    x_list = list()
    y_list = list()
    h_sat = [0., 0.]
    P_sat = 0.

    if T < support_point.RPHandler.TC:
        support_point.set_variable("T", T)
        support_point.set_variable("Q", 0)

        h_sat[0] = support_point.get_variable("h")

        support_point.set_variable("T", T)
        support_point.set_variable("Q", 1)

        h_sat[0] = support_point.get_variable("h")
        P_sat = support_point.get_variable("P")

    support_point.set_variable("P", P_max)
    support_point.set_variable("T", T)

    x_list.append(support_point.get_variable("h"))
    y_list.append(P_max)

    for h in h_range:

        if h_sat[0] <= h <= h_sat[1]:

            x_list.append(h)
            y_list.append(P_sat)

        else:

            support_point.set_variable("h", h)
            P_res = support_point.get_variable("P")

            if 0 < P_res <= P_max:
                x_list.append(h)
                y_list.append(P_res)

        pbar_iso_t.update(1)

    if len(x_list) > 0:
        isolines["T"].update({

            '{:.0f} C'.format(T): {

                "x": x_list,
                "y": y_list

            }

        })

pbar_iso_t.close()


# %%-------------------------------------   GET ISO-RHO LINES                   -------------------------------------> #
rho_range = get_np_array(150, 1000, 18)  # in kg/m^3
pbar_iso_h = tqdm(desc="iso-h lines calculation", total=n_isolines_points * len(rho_range))

for rho in rho_range:

    x_list = list()
    y_list = list()
    h_sat = [0., 0.]
    P_sat = 0.

    support_point.set_variable("rho", rho)

    for h in h_range:

        support_point.set_variable("h", h)
        P_res = support_point.get_variable("P")

        if P_res > 0:

            x_list.append(h)
            y_list.append(P_res)

        pbar_iso_h.update(1)

    if len(x_list) > 0:

        isolines["rho"].update({

            '{:.0f} kg/m^3'.format(rho): {

                "x": x_list,
                "y": y_list

            }

        })

pbar_iso_h.close()


# %%-------------------------------------   CALCULATE STANDARD BHE              -------------------------------------> #
depth = 1750
real_bhe = SimplifiedBHE(

    input_thermo_point=PlantThermoPoint(["CarbonDioxide"], [1]),
    dz_well=depth, t_rocks=110, use_rk=True

)

real_bhe.points[0].set_variable("T", 35)
real_bhe.points[0].set_variable("rho", 700)

real_bhe.update()

real_bhe_h = list()
real_bhe_P = list()

for point in real_bhe.points:
    real_bhe_h.append(point.get_variable("h"))
    real_bhe_P.append(point.get_variable("P"))

print(real_bhe)


# %%-------------------------------------   COMPARE WITH DISCRETIZATION         -------------------------------------> #
n_points = 1000
depth_list = np.linspace(0, depth, n_points, endpoint=True)
dz = depth_list[1] - depth_list[0]

down_point = real_bhe.points[0].duplicate()
up_point = real_bhe.points[2].duplicate()

p_up_list = list()
h_up_list = list()
rho_up_list = list()
p_down_list = list()
h_down_list = list()
rho_down_list = list()

for depth in depth_list[1:]:

    p_down = down_point.get_variable("P")
    h_down = down_point.get_variable("h")
    rho_down = down_point.get_variable("rho")

    p_up = up_point.get_variable("P")
    h_up = up_point.get_variable("h")
    rho_up = up_point.get_variable("rho")

    p_up_list.append(p_up)
    h_up_list.append(h_up)
    rho_up_list.append(rho_up)
    p_down_list.append(p_down)
    h_down_list.append(h_down)
    rho_down_list.append(rho_down)

    down_point.set_variable("P", p_down + rho_down * 9.81 * dz / 1e6)
    down_point.set_variable("h", h_down + 9.81 * dz / 1e3)

    up_point.set_variable("P", p_up - rho_up * 9.81 * dz / 1e6)
    up_point.set_variable("h", h_up - 9.81 * dz / 1e3)

disc_p_list = p_down_list + p_up_list
disc_h_list = h_down_list + h_up_list
disc_rho_list = rho_down_list + rho_up_list


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
highlight_isolines = False
plot_well_behaviour = True
plot_disc_well = True

fig, (ax_1) = plt.subplots(1, 1, dpi=200)
fig.set_size_inches(10, 5)

ax_1.plot(sat_h, sat_p, linewidth=2, color="black")

lines = list()
isoline_styles = {

    "T": "--",
    "rho": "-."

}

for line_key in isolines.keys():

    for key in isolines[line_key].keys():

        color = "black"
        width = 0.4
        style = isoline_styles[line_key]

        if highlight_isolines:

            if " C" in key:

                T_iso = float(key.strip(" C"))
                if T_iso == real_bhe.points[0].get_variable("T") or T_iso == real_bhe.points[2].get_variable("T"):

                    color = "white"
                    width = 0.8
                    style = "-"

            if " kg/m^3" in key:

                rho_iso = float(key.strip(" kg/m^3"))
                if rho_iso == real_bhe.points[0].get_variable("rho"):

                    color = "white"
                    width = 0.8
                    style = "-"

        lines.append(ax_1.plot(

            isolines[line_key][key]["x"], isolines[line_key][key]["y"],
            linewidth=width, color=color, alpha=0.5,
            linestyle=style, label=key

        )[0])


if plot_well_behaviour:

    ax_1.plot(real_bhe_h, real_bhe_P, "o", linewidth=2, markerfacecolor='none', markeredgecolor="white")

    if plot_disc_well:

        ax_1.plot(disc_h_list, disc_p_list, linewidth=2, color="white")

    else:

        ax_1.plot(real_bhe_h, real_bhe_P, linewidth=2, color="white")

contour = ax_1.contourf(h_mesh, p_mesh, res, n_points * 2)

ax_1.set_yscale("log")
ax_1.set_xlim(left=200, right=450)
ax_1.set_ylim(bottom=4, top=50)

ax_1.set_xlabel(

    r'$h\ [kJ/kg]$',
    fontsize='large', loc='right',

)

ax_1.set_ylabel(

    r'$P\ [MPa]$',
    fontsize='large', loc='top',

)

cb = fig.colorbar(contour)
cb.set_label(

    r'$log_{10}({c_o}^2)\ [-]$',
    fontsize='large', loc='top',

)

ax_1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.show()
