# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.power_plants.HTHP.subclasses.GEO_heat_pump_py import (

    WaterHeatPumpThermo,
    DirectWaterHeatPumpThermo

)
from main_code.power_plants.HTHP.subclasses.CO2_heat_pump_py import CO2HeatPumpThermo
from main_code.support.other.matplolib_stiles import BASE_COLORS, format_title_from_key
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%------------   INITIALIZE CALCULATIONS                -----------------------------------------------------------> #
depth = 800     # [m]
t_rock = 90     # [°C]
p_steam = 1     # [MPa]
t_in_BHE = 20   # [°C]

tmp_co2 = PlantThermoPoint(["CarbonDioxide"], [1])
tmp_h2o = PlantThermoPoint(["water"], [1])

tmp_co2.set_variable("T", t_in_BHE)
tmp_co2.set_variable("Q", 0)
t_in_BHE_co2 = tmp_co2.get_variable("P") + 0.1

tmp_h2o.set_variable("T", t_in_BHE)
tmp_h2o.set_variable("Q", 0)
t_in_BHE_h20 = tmp_h2o.get_variable("P") + 0.1

HTTP_dict = {

    "1) Direct sCO2 Heat Pump": CO2HeatPumpThermo(

        P_steam=p_steam,
        BHE_depth=depth, T_rock=t_rock,
        T_in_BHE=t_in_BHE, P_in_BHE=0.2, T_ambient=(t_in_BHE - 5),

    ),

    "2) Direct Steam Generation": DirectWaterHeatPumpThermo(

        P_steam=p_steam,
        BHE_depth=depth, T_rock=t_rock,
        T_in_BHE=t_in_BHE, P_in_BHE=t_in_BHE_h20, T_ambient=(t_in_BHE - 5),

    ),

    "3) Indirect Water Heat Pump\nnPentane": WaterHeatPumpThermo(

        P_steam=p_steam,
        BHE_depth=depth, T_rock=t_rock,
        T_in_BHE=t_in_BHE, P_in_BHE=t_in_BHE_h20, T_ambient=(t_in_BHE - 5),
        HTHP_fluid=["pentane"], HTHP_fluid_comp=[1]

    ),

    "4) Indirect Water Heat Pump\nWater": WaterHeatPumpThermo(

        P_steam=p_steam,
        BHE_depth=depth, T_rock=t_rock,
        T_in_BHE=t_in_BHE, P_in_BHE=t_in_BHE_h20, T_ambient=(t_in_BHE - 5),
        HTHP_fluid=["water"], HTHP_fluid_comp=[1]

    )

}


# %%------------   CALCULATION                            -----------------------------------------------------------> #
n_points = 20
result_dict = dict()

t_perc_lim = [0.1, 0.9]
range_lim = [5, 35]

t_sg_perc_list = np.linspace(t_perc_lim[0], t_perc_lim[1], n_points)
t_sc_perc_list = np.linspace(t_perc_lim[0], t_perc_lim[1], n_points)
range_eva_list = np.linspace(range_lim[0], range_lim[1], n_points)

pbar = tqdm(desc="calculation", total=n_points * len(HTTP_dict.keys()))

# 1 - Direct sCO2 Heat Pump
key = "1) Direct sCO2 Heat Pump"
HTTP = HTTP_dict[key]

COP_list = list()
m_ratio_list = list()

for t_sg_perc in t_sg_perc_list:

    HTTP.T_SG_perc = t_sg_perc
    HTTP.calculate(calculate_thermo_only=True)

    COP_list.append(HTTP.COP)
    m_ratio_list.append(2.45*HTTP.m_dot_ratio)

    pbar.update(1)

result_dict.update({

    key: {

        "x": t_sg_perc_list,
        "COP": COP_list,
        "m_ratio": m_ratio_list,
        "x_label": r'$T_{{SG\ \%}}\ [-]$',
        "x_lim": t_perc_lim,

    }

})

# 2 - Direct Steam Generation
key = "2) Direct Steam Generation"
HTTP = HTTP_dict[key]

COP_list = list()
m_ratio_list = list()

for t_sc_perc in t_sc_perc_list:

    HTTP.T_sc_perc = t_sc_perc
    HTTP.calculate(calculate_thermo_only=True)

    COP_list.append(HTTP.COP)
    m_ratio_list.append(HTTP.m_dot_ratio)

    pbar.update(1)

result_dict.update({

    key: {

        "x": t_sc_perc_list,
        "COP": COP_list,
        "m_ratio": m_ratio_list,
        "x_label": r'$T_{{SG\ \%}}\ [-]$',
        "x_lim": t_perc_lim,

    }

})

# 3 - Indirect Water Heat Pump - nPentane
key = "3) Indirect Water Heat Pump\nnPentane"
HTTP = HTTP_dict[key]

COP_list = list()
m_ratio_list = list()

for range_eva in range_eva_list:

    HTTP.range_eva = range_eva
    HTTP.calculate(calculate_thermo_only=True)

    COP_list.append(HTTP.COP)
    m_ratio_list.append(HTTP.m_dot_ratio)

    pbar.update(1)

result_dict.update({

    key: {

        "x": range_eva_list,
        "COP": COP_list,
        "m_ratio": m_ratio_list,
        "x_label": r'$range_{{EVA}}\ [°C]$',
        "x_lim": range_lim,

    }

})

# 4 - Indirect Water Heat Pump - Water
key = "4) Indirect Water Heat Pump\nWater"
HTTP = HTTP_dict[key]
range_eva_list = np.linspace(5, 35, 20)

COP_list = list()
m_ratio_list = list()

for range_eva in range_eva_list:

    HTTP.range_eva = range_eva
    HTTP.calculate(calculate_thermo_only=True)

    COP_list.append(HTTP.COP)
    m_ratio_list.append(HTTP.m_dot_ratio)

    pbar.update(1)

result_dict.update({

    key: {

        "x": range_eva_list,
        "COP": COP_list,
        "m_ratio": m_ratio_list,
        "x_label": r'$range_{{EVA}}\ [°C]$',
        "x_lim": range_lim,

    }

})

pbar.close()


# %%------------   DATA PLOT                              -----------------------------------------------------------> #
# Figure initialization
fig, ax_list = plt.subplots(2, 2, dpi=150)
fig.set_size_inches(12, 8)

# Math Test Font Definition
plt.rcdefaults()


i = 0
for key in HTTP_dict.keys():

    main_ax = ax_list[np.floor_divide(i, 2)][np.mod(i, 2)]
    sec_ax = main_ax.twinx()

    x_list = result_dict[key]["x"]
    y_list = result_dict[key]["COP"]
    y_2_list = result_dict[key]["m_ratio"]

    line1 = main_ax.plot(

        x_list, y_list,
        label='COP',
        color=BASE_COLORS[0]

    )

    line2 = sec_ax.plot(

        x_list, y_2_list,
        label=r'$m_{{ratio}}$',
        color=BASE_COLORS[1]

    )

    # add legend
    if i == 0:

        lns = line1 + line2
        labs = [l.get_label() for l in lns]

        main_ax.legend(lns, labs, loc='upper left')

    i = i + 1

    # Axes Title
    main_ax.set_xlabel(

        result_dict[key]["x_label"],
        fontsize='large', loc='right',
        weight='light'

    )
    main_ax.set_ylabel(

        "COP [-]",
        fontsize='large', loc='top',
        weight='light'

    )
    sec_ax.set_ylabel(

        r'$m_{{ratio}}\ [-]$',
        fontsize='large', loc='top',
        weight='light'

    )

    # Grid Style
    main_ax.tick_params(axis=u'both', which=u'minor', length=0)
    main_ax.grid(which='major', linewidth=0.8, alpha=0.75)

    # Set x axis lim
    main_ax.set_xlim(result_dict[key]["x_lim"])
    sec_ax.set_xlim(result_dict[key]["x_lim"])

    # Set Title
    main_ax.set_title(format_title_from_key(key), loc="left", fontsize=10, fontfamily="sans-serif")

plt.tight_layout(pad=2)
plt.show()


# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #
current_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "4 - Application Case Studies", "4 - HTHP"

)
output_directory = os.path.join(current_folder, "outputs")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "COP and m_ratio.png"))

