# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.support_functions import get_np_array
from main_code.support.other.label_lines import label_lines
import matplotlib.pyplot as plt
from tqdm import tqdm


# %%-------------------------------------   INITIALIZATION                      -------------------------------------> #
CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", 45)
CO2_input.set_variable("rho", 800)

water_input = PlantThermoPoint(["Water"], [1])
water_input.set_variable("T", 45)
water_input.set_variable("P", 0.2)

bhe_CO2 = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=1000, t_rocks=150

)

bhe_water = SimplifiedBHE(

    input_thermo_point=water_input,
    dz_well=1000, t_rocks=150

)

n_points = 100
depth_points = get_np_array(500, 5000, n_points)    # in m
t_in_points = get_np_array(90, 230, 8)              # in °C
res_dict = dict()

str_form = r'$T_{{rock}}={:.0f}°C$'
line_styles = {

    str_form.format(90): {

        "style": "-",
        "size": 2

    },

    str_form.format(110): {

        "style": "--",
        "size": 1

    },

    str_form.format(130): {

        "style": "--",
        "size": 1.5

    },

    str_form.format(150): {

        "style": "--",
        "size": 2

    },

    str_form.format(170): {

        "style": "-.",
        "size": 1

    },

    str_form.format(190): {

        "style": "-.",
        "size": 1.5

    },

    str_form.format(210): {

        "style": "-.",
        "size": 2

    },

    str_form.format(230): {

        "style": "-",
        "size": 1

    },

    "default": {

        "style": "-",
        "size": 0.5

    },

}
line_styles.setdefault("default")


# %%-------------------------------------   CO2 CALCULATION                     -------------------------------------> #
pbar = tqdm(desc="evaluating CO2 points", total=n_points * len(t_in_points))

for i in range(len(t_in_points)):

    x_CO2 = list()
    y_CO2 = list()

    for j in range(len(depth_points)):

        bhe_CO2.dz_well = depth_points[j]
        bhe_CO2.t_rocks = t_in_points[i]

        bhe_CO2.update()

        if 0 <= bhe_CO2.eta_II <= 1:

            x_CO2.append(depth_points[j])
            y_CO2.append(bhe_CO2.eta_II)

            pbar.update(1)

        else:

            pbar.update(len(depth_points) - j)
            break

    key = str_form.format(t_in_points[i])

    if key not in res_dict.keys():
        res_dict.update({

            key: {

                "Water": dict(),
                "CO2": dict()

            }
        })

    res_dict[key]["CO2"].update({

        "x": x_CO2,
        "y": y_CO2

    })

pbar.close()


# %%-------------------------------------   WATER CALCULATION                   -------------------------------------> #
pbar = tqdm(desc="evaluating water points", total=n_points * len(t_in_points))

for i in range(len(t_in_points)):

    x_water = list()
    y_water = list()

    for j in range(len(depth_points)):

        bhe_water.dz_well = depth_points[j]
        bhe_water.t_rocks = t_in_points[i]
        bhe_water.update()

        if 0 <= bhe_water.eta_II <= 1:

            x_water.append(depth_points[j])
            y_water.append(bhe_water.eta_II)
            pbar.update(1)

        else:

            pbar.update(len(depth_points) - j)
            break

    key = str_form.format(t_in_points[i])

    if key not in res_dict.keys():

        res_dict.update({

            key: {

                "Water": dict(),
                "CO2": dict()

            }
        })

    res_dict[key]["Water"].update({

        "x": x_water,
        "y": y_water

    })

pbar.close()


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
fig, (ax) = plt.subplots(1, 1, dpi=200)
fig.set_size_inches(10, 5)

lines = list()

label_list = list(res_dict.keys())

# Print CO2
for label in res_dict.keys():

    lines.append(ax.plot(

        res_dict[label]["CO2"]["x"], res_dict[label]["CO2"]["y"],
        color="black", label=label,
        linestyle=line_styles[label]["style"],
        linewidth=line_styles[label]["size"]

    )[0])

# Print Water
for label in [label_list[0], label_list[-1]]:

    lines.append(ax.plot(

        res_dict[label]["Water"]["x"], res_dict[label]["Water"]["y"],
        color='#86b1f7', label=label,
        linewidth=2

    )[0])

# Axes Title
ax.set_xlabel(

    r'$Depth_{{\ HE\ Section}}\ [m]$',
    fontsize='large', loc='right',
    weight='light'

)

ax.set_ylabel(

    r'$Exergy\ Efficiency\ [-]$',
    fontsize='large', loc='top',
    weight='light'

)

ax.set_xlim(left=500, right=5000)
ax.set_ylim(bottom=0.4, top=0.75)

x_lab_pos = [

    1850,
    2400,
    2950,
    3500,
    4150,
    4750,
    3900,
    4750,

]
for i in range(2):
    x_lab_pos.append(4750)

# lines[-3].set_label("_line")
label_lines(lines, xvals=x_lab_pos, fontsize=8)

ax.grid(which='major', linewidth=0.8)
plt.show()
