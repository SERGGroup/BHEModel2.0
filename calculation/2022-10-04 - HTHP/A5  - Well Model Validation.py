# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.simplified_well.heating_sections.subclasses.EGS_heating_section import EGSHeatingSection
from main_code.well_model.simplified_well.heating_sections.subclasses.default_class import DefaultHeatingSection
from main_code.well_model.geometry_based_well_models.simple_pressure_losses_model import PressureLossesBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.matplolib_stiles import ColorFader
from openpyxl import load_workbook
from main_code import constants
import matplotlib.pyplot as plt
import numpy as np
import os


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
    T_rock = t_amb + t_grad * dz_well / 1e3

    if pressure_losses:

        d_inj = 0.41
        d_prod = 0.27

    else:

        d_inj = None
        d_prod = None

    bhe_in = PressureLossesBHE(

        input_thermo_point=tmp_co2_input,
        dz_well=dz_well, t_rocks=T_rock, use_rk=True,
        calc_discrete_pressure_losses=discrete_losses,
        d_inj=d_inj, d_prod=d_prod

    )

    if res_losses:

        EGSHeatingSection(bhe_in, R_res=R_res, Q_res=Q_res)

    else:

        DefaultHeatingSection(bhe_in)

    bhe_in.update()

    if not pump_power == 0.:

        w_iso = pump_power / eta_pump
        h_cond_iso = tmp_co2_input.get_variable("h") - w_iso / tmp_co2_input.m_dot
        s_iso = tmp_co2_input.get_variable("s")

        tmp_co2_input.set_variable("s", s_iso)
        tmp_co2_input.set_variable("h", h_cond_iso)

    p_cond_out = tmp_co2_input.get_variable("P")
    turb_output = bhe_in.points[-1].duplicate()
    turb_output.set_to_expansion_result(p_cond_out, eta_turb, bhe_in.points[-1])

    dp_turb = bhe_in.points[-1].get_variable("P") - p_cond_out
    dh_turb = bhe_in.points[-1].get_variable("h") - turb_output.get_variable("h")
    W_turb = dh_turb * turb_output.m_dot
    Q_res = bhe_in.Q_bottom

    return dp_turb, dh_turb, W_turb, Q_res, bhe_in


# %%------------   CHECK ADAMS FROM EXCEL                 -----------------------------------------------------------> #
file_directory = os.path.join(constants.CALCULATION_FOLDER, "2022-11-21 - Model Validation", "res")
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

        try:

            sheet[row][col_w_turb+j].value = W_turb
            sheet[row][col_dp_turb+j].value = dp_turb

        except:

            sheet[row][col_w_turb+j].value = 0.
            sheet[row][col_dp_turb+j].value = 0.

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


# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #
current_folder = os.path.join(os.path.dirname(constants.RES_FOLDER), "2022-10-04 - HTHP")
output_directory = os.path.join(current_folder, "outputs")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "validation.png"))
