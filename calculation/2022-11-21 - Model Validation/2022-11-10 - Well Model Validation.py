# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.simplified_well.heating_sections.subclasses.EGS_heating_section import EGSHeatingSection
from main_code.simplified_well.heating_sections.subclasses.default_class import DefaultHeatingSection
from main_code.simplified_well.simplified_well_subclasses import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.matplolib_stiles import ColorFader
from matplotlib.ticker import MultipleLocator
from openpyxl import load_workbook
from main_code import constants
import matplotlib.pyplot as plt
from matplotlib import ticker
from tqdm import tqdm
import numpy as np
import os


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

eta_pump = 0.9
eta_turb = 0.78

CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("Q", 0)

p_sat = CO2_input.get_variable("P")
CO2_input.set_variable("T", T_amb + dT_appr)
CO2_input.set_variable("P", p_sat + dp_sat / 1e3)


def evaluate_turbine_power(

        co2_in, t_amb, t_grad, dz_well,
        pressure_losses=False, res_losses=False,
        discrete_losses=False, pump_power=0.,
        R_res=None, Q_res=None

):

    tmp_co2_input = co2_in.duplicate()

    if not pump_power == 0.:

        h_out = tmp_co2_input.get_variable("h") + pump_power / tmp_co2_input.m_dot

        w_iso = pump_power / eta_pump
        h_iso = tmp_co2_input.get_variable("h") + w_iso / tmp_co2_input.m_dot
        s_iso = tmp_co2_input.get_variable("s")

        tmp_co2_input.set_variable("s", s_iso)
        tmp_co2_input.set_variable("h", h_iso)

        p_out = tmp_co2_input.get_variable("p")

        tmp_co2_input.set_variable("p", p_out)
        tmp_co2_input.set_variable("h", h_out)

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

        EGSHeatingSection(bhe_in, R_res=R_res, Q_res=Q_res)

    else:

        DefaultHeatingSection(bhe_in)

    bhe_in.update()

    turb_output = bhe_in.points[-1].duplicate()
    turb_output.set_to_expansion_result(p_sat + dp_sat / 1e3, eta_turb, bhe_in.points[-1])

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

T_grad = 50         # [°C/km]
depth = 1000        # [m]

R_res = 22.5
q_res = 8263

m_dot = 49             # [kg/s]
pump_power = 0.0       # [kWe]

CO2_input.m_dot = m_dot

dp_turb, dh_turb, W_turb, q_res, bhe_out = evaluate_turbine_power(

    CO2_input, T_amb, T_grad, depth,
    pump_power=pump_power, res_losses=False,
    pressure_losses=True, R_res=R_res, Q_res=q_res

)

# %%------------   CHECK ADAMS FROM EXCEL                 -----------------------------------------------------------> #

file_directory = os.path.join(os.path.dirname(constants.RES_FOLDER), "2022-11-21 - Model Validation", "res")
file_name = os.path.join(file_directory, "1-s2.0-S0306261914012124-mmc2.xlsx")
wb = load_workbook(filename=file_name)
sheet = wb.worksheets[0]
results = [list(), list(), list(), list()]

for row in range(7, 17):

    T_grad = abs(sheet[row][3].value)
    depth = abs(sheet[row][4].value)
    R_res = abs(sheet[row][5].value)
    q_res = abs(sheet[row][6].value)
    m_dot = abs(sheet[row][7].value)
    pump_power = abs(sheet[row][8].value)

    CO2_input.m_dot = m_dot

    j = 0
    col_w_turb = 11
    col_dp_turb = 15

    results[0].append(sheet[row][10].value)

    for losses in [[True, True], [True, False], [False, False]]:

        dp_turb, dh_turb, W_turb, q_res_calc, bhe_out = evaluate_turbine_power(

            CO2_input, T_amb, T_grad, depth,
            pump_power=pump_power, res_losses=losses[0],
            pressure_losses=losses[1], R_res=R_res, Q_res=q_res

        )

        sheet[row][col_w_turb+j].value = W_turb
        sheet[row][col_dp_turb+j].value = dp_turb
        results[j + 1].append(sheet[row][col_w_turb+j].value)
        j += 1

wb.save(file_name)

results[0].sort()
results[1].sort()

results[0] = np.array(results[0])*1.08


# %%------------   PLOT ADAMS EXCEL RESULTS               -----------------------------------------------------------> #

x_label = "Wturb (Adams et al.) [kWe]"
y_label = "Wturb (current study) [kWe]"
cf = ColorFader()

fig, (ax_1) = plt.subplots(1, 1, dpi=200)
fig.set_size_inches(10, 5)

i = 0
linestyles = ['solid', 'dashed', 'dotted', 'dashdot']

ax_1.set_xscale("log")
ax_1.set_yscale("log")
ax_1.grid(visible=True, which="major")
ax_1.grid(visible=True, which="minor", alpha=0.2)

ax_1.plot(

    results[0], results[0], linestyle=linestyles[3],
    color=cf.get_color(0), alpha=0.5, linewidth=1

)

error = 0.15
ax_1.fill_between(

    results[0], np.array(results[0])*(1+error), y2=np.array(results[0])*(1-error),
    color=cf.get_color(0), alpha=0.3

)

ax_1.plot(

    results[0], results[1], linestyle=linestyles[0],
    color=cf.get_color(1), linewidth=2

)

ax_1.set_xlabel(

    x_label,
    fontsize='large', loc='right',

)

ax_1.set_ylabel(

    y_label,
    fontsize='large', loc='top',

)

plt.show()


# %%------------   CHECK ADAMS DICTIONARY                 -----------------------------------------------------------> #

""" 
    DATA FROM:

    Adamns et Al.
    "A comparison of electric power output of CO2 Plume Geothermal (CPG) and brine
     geothermal systems for varying reservoir conditions, Applied Energy"

"""

depth = 2500  # [m]
T_grad = 35  # [°C/km]
m_dot_list = np.linspace(40, 180, 100)
results = dict()

keys = ["no losses", "0 - Resources losses", "all losses", "discrete losses"]
pbar = tqdm(desc="profile calculation", total=len(m_dot_list) * len(keys))

for key in keys:

    discrete_losses = False
    pressure_losses = False
    res_losses = False

    if key == "0 - Resources losses":
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
        dp_turb, dh_turb, W_turb, q_res, bhe_out = evaluate_turbine_power(

            CO2_input, T_amb, T_grad, depth,
            pressure_losses=pressure_losses,
            res_losses=res_losses

        )

        dp_list.append(dp_turb)
        dh_list.append(dh_turb)
        W_turb_list.append(W_turb)
        Q_res_list.append(q_res)
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

# %%------------   PLOT ADAMS DICTIONARY RESULTS          -----------------------------------------------------------> #

x_list = "m_dot"
x_label = r'$m_{dot}\ [kg/s]$'

y_list = "dp_list"
y_label = r'$dP_{turb}\ [MPa]$'

# y_list = "dh_list"
# y_label = r'$dh_{turb}\ [kJ/kg]$'

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
ax_1.set_xlim(left=40, right=180)


ax_1.xaxis.set_major_locator(MultipleLocator(20))
ax_1.xaxis.set_minor_locator(MultipleLocator(5))
ax_1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

ax_1.yaxis.set_minor_locator(MultipleLocator(1))
ax_1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax_1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

ax_1.spines["right"].set_visible(False)
ax_1.spines["top"].set_visible(False)

plt.show()
