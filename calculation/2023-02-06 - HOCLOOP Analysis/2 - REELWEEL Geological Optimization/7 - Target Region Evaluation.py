# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.gradient_fill import gradient_fill
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from main_code import constants
import numpy as np
import os

# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

h_try = 375
p_list = np.linspace(10, 50, 50)

t_pure_list = list()
t_mix_list = list()

t_h_pure_list = list()
t_h_mix_list = list()

tmp_point = PlantThermoPoint(["Carbon Dioxide"], [1])
tmp_point_mix = PlantThermoPoint(["Carbon Dioxide", "Nitrogen"], [0.9, 0.1])

for p in p_list:

    def opt_function(h):

        tmp_point.set_variable("P", p)
        tmp_point.set_variable("h", h)

        return -tmp_point.get_variable("cp")

    h = minimize_scalar(opt_function, method='Bounded', bounds=(250, 450)).x
    print(h)

    tmp_point.set_variable("P", p)
    tmp_point.set_variable("h", h)
    t_pure_list.append(tmp_point.get_variable("T"))

    tmp_point.set_variable("P", p)
    tmp_point.set_variable("h", h_try)
    t_h_pure_list.append(tmp_point.get_variable("T"))


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #

for p in p_list:
    def opt_function(h):

        tmp_point_mix.set_variable("P", p)
        tmp_point_mix.set_variable("h", h)

        return -tmp_point_mix.get_variable("cp")

    h = minimize_scalar(opt_function, method='Bounded', bounds=(250, 450)).x
    print(h)

    tmp_point_mix.set_variable("P", p)
    tmp_point_mix.set_variable("h", h)
    t_mix_list.append(tmp_point_mix.get_variable("T"))

    tmp_point_mix.set_variable("P", p)
    tmp_point_mix.set_variable("h", h_try)
    t_h_mix_list.append(tmp_point_mix.get_variable("T"))


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #

fig = plt.figure(dpi=150)
fig.set_size_inches(10, 8)
ax = fig.add_subplot(1, 1, 1)


line_mix, _ = gradient_fill(

    np.array(p_list), np.array(t_mix_list),
    ax=ax, fill_dir="top", grad_dir="S",
    color="gold"

)

line_pure, _ = gradient_fill(

    np.array(p_list), np.array(t_pure_list),
    ax=ax, fill_dir="top", grad_dir="S",
    color="goldenrod"

)

ax.plot(p_list, t_h_mix_list, "--", color=line_mix.get_color(), zorder=100)
line_dashed = ax.plot(p_list, t_h_pure_list, "--", color="silver")
ax.plot(p_list, t_h_pure_list, "--", color=line_pure.get_color(), zorder=100)

ax.set_ylabel(r'$T_{heel}$ [°C]', fontsize='large', loc='center')
ax.set_xlabel(r'$P_{heel}$ [MPa]', fontsize='large', loc='center')
ax.set_xlabel(r'$P_{heel}$ [MPa]', fontsize='large', loc='center')
ax.set_title("Pseudo-Boiling Line", fontsize='xx-large', loc='center')

ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
plt.legend([line_pure, line_mix, line_dashed[0]], ["Pure $CO_2$", "90% $CO_2$, 10% $N_2$", "h = 375 [kJ/kg]"])
plt.show()

# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #

RES_FOLDER = os.path.join(constants.CALCULATION_FOLDER, "2023-02-06 - HOCLOOP Analysis", "0 - Resources")
output_directory = os.path.join(RES_FOLDER, "output")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "7 - Target Region Evaluation.png"))