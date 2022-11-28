# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.simplified_BHE.heating_sections.subclasses.EGS_heating_section import EGSHeatingSection
from main_code.simplified_BHE.heating_sections.subclasses.default_class import DefaultHeatingSection
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.simplified_BHE.simplified_BHE import SimplifiedBHE
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from matplotlib import ticker
from tqdm import tqdm
import numpy as np


# %%------------   CHECK PREUSS                           -----------------------------------------------------------> #

""" 
    DATA FROM:

    Pruess, K.
    "Enhanced geothermal systems (EGS) using CO2 as working fluid
     A novel approach for generating renewable energy with simultaneous sequestration of carbon"

"""

P_in = 5.74  # [MPa]
T_in = 20  # [°C]

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_in)
CO2_input.set_variable("P", P_in)

water_input = PlantThermoPoint(["Water"], [1])
water_input.set_variable("T", T_in)
water_input.set_variable("P", P_in)

bhe_CO2 = SimplifiedBHE(

    input_thermo_point=CO2_input,
    dz_well=5000, T_rocks=200

)

bhe_water = SimplifiedBHE(

    input_thermo_point=water_input,
    dz_well=5000, T_rocks=200

)

bhe_CO2.update()
bhe_water.update()

P_CO2_down = bhe_CO2.points[1].get_variable("P")
P_water_down = bhe_water.points[1].get_variable("P")

P_CO2_up = bhe_CO2.points[3].get_variable("P")
P_water_up = bhe_water.points[3].get_variable("P")

# %%------------   CHECK ADAMS INIT                       -----------------------------------------------------------> #

T_amb = 15  # [°C]
dT_appr = 7  # [°C]
dp_sat = 50  # [kPa]

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("Q", 0)

p_sat = CO2_input.get_variable("P")
CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("P", p_sat + dp_sat / 1e3)


def evaluate_turbine_power(

        co2_in, t_amb, t_grad, dz_well,
        pressure_losses=False, res_losses=False,
        discrete_losses=False

):

    tmp_co2_input = co2_in.duplicate()
    T_rock = t_amb + t_grad * dz_well / 1e3

    if pressure_losses:

        d_inj = 0.41
        d_prod = 0.27

    else:

        d_inj = None
        d_prod = None

    bhe_in = SimplifiedBHE(

        input_thermo_point=tmp_co2_input,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True,
        discretize_p_losses=discrete_losses,
        d_inj=d_inj, d_prod=d_prod

    )

    if res_losses:

        EGSHeatingSection(bhe_in)

    else:

        DefaultHeatingSection(bhe_in)

    bhe_in.update()

    turb_output = bhe_in.points[-1].duplicate()
    turb_output.set_to_expansion_result(p_sat + dp_sat / 1e3, 0.78, bhe_in.points[-1])

    dp_turb = bhe_in.points[-1].get_variable("P") - (p_sat + dp_sat / 1e3)
    dh_turb = bhe_in.points[-1].get_variable("h") - turb_output.get_variable("h")
    W_turb = dh_turb * turb_output.m_dot
    Q_res = bhe_in.Q_bottom

    return dp_turb, dh_turb, W_turb, Q_res, bhe_in


# %%------------   CHECK ADAMS SINGLE                     -----------------------------------------------------------> #

""" 
    DATA FROM:

    Adamns et Al.
    "A comparison of electric power output of CO2 Plume Geothermal (CPG) and brine
     geothermal systems for varying reservoir conditions, Applied Energy"

"""

m_dot = 88  # [kg/s]
depth = 1500  # [m]
T_grad = 35  # [°C/km]

dp_turb, dh_turb, W_turb, Q_res, bhe_out = evaluate_turbine_power(CO2_input, T_amb, T_grad, depth)

# %%------------   CHECK ADAMS DICTIONARY                 -----------------------------------------------------------> #

""" 
    DATA FROM:

    Adamns et Al.
    "A comparison of electric power output of CO2 Plume Geothermal (CPG) and brine
     geothermal systems for varying reservoir conditions, Applied Energy"

"""

depth = 2500  # [m]
T_grad = 35  # [°C/km]
m_dot_list = np.linspace(1, 200, 100)
results = dict()

keys = ["no losses", "res losses", "all losses", "discrete losses"]
pbar = tqdm(desc="profile calculation", total=len(m_dot_list) * len(keys))

for key in keys:

    discrete_losses = False
    pressure_losses = False
    res_losses = False

    if key == "res losses":
        res_losses = True

    if key == "all losses":
        res_losses = True
        pressure_losses = True

    if key == "discrete losses":
        res_losses = True
        pressure_losses = True
        discrete_losses = True

    dp_list = list()
    dh_list = list()
    W_turb_list = list()
    Q_res_list = list()
    bhe_list = list()

    for m_dot in m_dot_list:

        CO2_input.m_dot = m_dot
        dp_turb, dh_turb, W_turb, Q_res, bhe_out = evaluate_turbine_power(

            CO2_input, T_amb, T_grad, depth,
            pressure_losses=pressure_losses,
            res_losses=res_losses

        )

        dp_list.append(dp_turb)
        dh_list.append(dh_turb)
        W_turb_list.append(W_turb)
        Q_res_list.append(Q_res)
        bhe_list.append(bhe_out)

        pbar.update(1)

    results.update({

        key: {

            "dp_list": dp_list,
            "dh_list": dh_list,
            "W_turb_list": W_turb_list,
            "Q_res_list": Q_res_list,
            "bhe_list": bhe_list

        }

    })

pbar.close()

# %%------------   PLOT ADAMS CHECK RESULTS               -----------------------------------------------------------> #

x_list = "m_dot"
x_label = r'$m_{dot}\ [kg/s]$'

# y_list = "dp_list"
# y_label = r'$dP_{turb}\ [MPa]$'

y_list = "dh_list"
y_label = r'$dh_{turb}\ [kJ/kg]$'

fig, (ax_1) = plt.subplots(1, 1, dpi=200)
fig.set_size_inches(10, 5)

i = 0
linestyles = ['dashed', 'dotted', 'dashdot', 'solid']

for key in results.keys():

    sub_dict = results[key]

    if x_list == "m_dot":

        x_array = m_dot_list

    else:

        x_array = sub_dict[x_list]

    ax_1.plot(x_array, sub_dict[y_list], linestyle=linestyles[i], color='black')
    i += 1

ax_1.set_xlabel(

    x_label,
    fontsize='large', loc='right',

)

ax_1.set_ylabel(

    y_label,
    fontsize='large', loc='top',

)

ax_1.set_ylim(bottom=0)
ax_1.set_xlim(left=0, right=200)


ax_1.xaxis.set_major_locator(MultipleLocator(25))
ax_1.xaxis.set_minor_locator(MultipleLocator(5))
ax_1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

ax_1.yaxis.set_minor_locator(MultipleLocator(1))
ax_1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax_1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

ax_1.spines["right"].set_visible(False)
ax_1.spines["top"].set_visible(False)

plt.show()
