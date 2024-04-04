# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from EESConnect import EESConnector
from main_code import constants
from datetime import datetime
from pandas import DataFrame
import numpy as np
import os.path


# %%------------   CALCULATIONS DATA                      -----------------------------------------------------------> #
dt_he = 40
t_turb_in = 20
eta_comp = 0.8
eta_turb = 0.75

n_T_max = 10
n_grad_T = 6
n_m_dot = 11

t_max_arr = np.linspace(80, 120, n_T_max)
grad_T_arr = np.linspace(50, 75, n_grad_T)
m_dot_well_arr = np.linspace(5, 15, n_m_dot)

t_max_arr, grad_T_arr, m_dot_well_arr = np.meshgrid(t_max_arr, grad_T_arr, m_dot_well_arr, indexing='ij')
inputs_names = ["t_max", "dt_he", "t_turb_in", "eta_comp", "eta_turb", "grad_T", "m_dot_well"]

calculation_dict = {}

for i in range(n_T_max):

    for j in range(n_grad_T):

        for k in range(n_m_dot):

            calculation_dict.update({

                "{}-{}-{}".format(i, j, k):  [

                    t_max_arr[i, j, k], dt_he, t_turb_in,
                    eta_comp, eta_turb, grad_T_arr[i, j, k], m_dot_well_arr[i, j, k]

                ]
            })


# %%------------   CALCULATE                              -----------------------------------------------------------> #
RES_FOLDER = os.path.join(

    constants.CALCULATION_FOLDER, "HOCLOOP Calculation",
    "WP4", "3 - ECOS 2024 Paper", "2 - EES calculation",

)

EES_FILE = os.path.join(RES_FOLDER, "0 - EES Files", "base heat pump - python.EES")

with EESConnector(EES_FILE, ees_decimal_separator=",", display_progress_bar=True) as ees:

    try:
        result = ees.calculate(calculation_dict)

    except:
        result = {}

print("Done!")

now = datetime.now()
now_str = "{}-{:02}-{:02} - {:02}.{:02}".format(now.year, now.month, now.day, now.hour, now.minute)
file_path = os.path.join(RES_FOLDER, "00 - Results", "calculation results files", "{}.xlsx".format(now_str))
df_results = list()

results_names = ["LCOH", "Q_DH", "w_net", "t_max_reach"]

for i in range(n_T_max):

    for j in range(n_grad_T):

        for k in range(n_m_dot):

            sub_list = calculation_dict["{}-{}-{}".format(i, j, k)]
            sub_list.extend(result["{}-{}-{}".format(i, j, k)])

            df_results.append(sub_list)

names = inputs_names
names.extend(results_names)
df = DataFrame(df_results)
df.columns = names
df.to_excel(file_path, index=False, startrow=1, startcol=1)
