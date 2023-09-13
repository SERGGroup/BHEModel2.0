# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.ground_models.FreeFemAnalysis import (

    FreeFEMAnalyzer, TimeDependentOptions, MeshOptions, ProblemDefinitionOptions

)
from main_code.support.other.excel_exporter import write_excel_sheet
from main_code.constants import CALCULATION_FOLDER, os
from collections.abc import Iterable
import matplotlib.pyplot as plt
import numpy as np


# %%------------   DEFINE WORKING FOLDERS                 -----------------------------------------------------------> #
calculation_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "1 - mesh sensitivity",
    "res", "calculation folder"

)

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "1 - mesh sensitivity",
    "res", "results"

)


# %%------------   CALCULATIONS                         -----------------------------------------------------------> #
a = 4
results = dict()
DT = 40.
for b in [5, 10, 15]:

    for multiplier in [1, 2, 3, 4]:

        print("{}x - {}%".format(multiplier, round(100 / b)))

        mo = MeshOptions(

            n_points=a*multiplier, n_points_circle=a*b*multiplier, graph_r_ratio=15,
            mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
            retrieve_mesh=False, add_visualization_mesh=True

        )

        pdo = ProblemDefinitionOptions(

            grad_rock=0.01, DT=DT,
            time_range=[1e-3, 1000], time_steps=60,
            n_plot=0, n_save=50

        )

        tdo = TimeDependentOptions(mesh_options=mo, problem_definition_options=pdo, mesh_only=False)
        ffa = FreeFEMAnalyzer(options=tdo, calculation_folder=calculation_folder)
        gross_result = ffa.calculate()

        time_list = list()
        value_list = list()

        for line in gross_result:

            line_split = line.split(";")
            time_list.append(float(line_split[0].strip()))
            value_list.append(float(line_split[1].strip()))

        results.update({

            "{}x - {}%".format(multiplier, round(100 / b)): {

                "time": time_list,
                "value": value_list,
                "calc_time": ffa.mean_calculation_time

            }

        })

        print("{}x - {}%".format(multiplier, round(100 / b)))


# %%------------   EVALUATE RELATIVE DIFFERENCES          -----------------------------------------------------------> #
res_keys_list = list(results.keys())
difference_dataframe = {

    'Mesh': {"unit": [None], "values": [res_keys_list[1:]]},
    'relative difference': {"unit": list(), "values": list()},
    'final difference': {"unit": list(), "values": list()},

}

log_difference_dataframe = {

    'Mesh': {"unit": [None], "values": [res_keys_list[1:]]},
    'relative difference': {"unit": list(), "values": list()},
    'final difference': {"unit": list(), "values": list()},

}

for i in range(len(res_keys_list)):

    key_i = res_keys_list[i]
    arr_i = np.array(results[key_i]["value"])
    log_arr_i = np.log10(arr_i)

    difference_dataframe['relative difference']["unit"].append(key_i)
    difference_dataframe['relative difference']["values"].append(list())
    difference_dataframe['final difference']["unit"].append(key_i)
    difference_dataframe['final difference']["values"].append(list())

    log_difference_dataframe['relative difference']["unit"].append(key_i)
    log_difference_dataframe['relative difference']["values"].append(list())
    log_difference_dataframe['final difference']["unit"].append(key_i)
    log_difference_dataframe['final difference']["values"].append(list())

    for j in range(i+1, len(res_keys_list)):

        key_j = res_keys_list[j]
        arr_j = np.array(results[key_j]["value"])
        log_arr_j = np.log10(arr_j)

        # Standard Differences
        d_arr_sqr = (arr_j - arr_i) ** 2
        rel_diff = np.sqrt(np.mean(d_arr_sqr / arr_i))
        d_final = np.abs((arr_j[-1] - arr_i[-1])/arr_i[-1])

        difference_dataframe['relative difference']["values"][-1].append(rel_diff)
        difference_dataframe['final difference']["values"][-1].append(d_final)


        # Log Based Differences
        d_arr_sqr = (log_arr_j - log_arr_i) ** 2
        rel_diff = np.sqrt(np.mean(d_arr_sqr / log_arr_i))
        d_final = np.abs((log_arr_j[-1] - log_arr_i[-1]) / log_arr_i[-1])

        log_difference_dataframe['relative difference']["values"][-1].append(rel_diff)
        log_difference_dataframe['final difference']["values"][-1].append(d_final)


# %%------------   SAVE RESULTS                           -----------------------------------------------------------> #
file_path = os.path.join(result_folder, "mesh-sensitivity-{}.xlsx".format(a))
timing_dataframe = {

    'Mesh': {"unit": [None], "values": [list()]},
    'Computation Time': {"unit": ["ms"], "values": [list()]},

}
for key in results.keys():

    timing_dataframe['Mesh']["values"][0].append(key)
    timing_dataframe['Computation Time']["values"][0].append(results[key]["calc_time"])

    dataframe = {

        'Time': {"unit": ["-"], "values": [results[key]["time"]]},
        'value': {"unit": ["-"], "values": [results[key]["value"]]},

    }

    write_excel_sheet(excel_path=file_path, sheet_name=key, data_frame=dataframe, overwrite="hard")

write_excel_sheet(excel_path=file_path, sheet_name="calc_time", data_frame=timing_dataframe, overwrite="hard")
write_excel_sheet(excel_path=file_path, sheet_name="diff", data_frame=difference_dataframe, overwrite="hard")
write_excel_sheet(excel_path=file_path, sheet_name="log_diff", data_frame=log_difference_dataframe, overwrite="hard")


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
def calculate_q(td):

    if issubclass(type(td), Iterable):

        result = list()
        for td_item in td:
            result.append(calculate_q(td_item))

        return result

    else:

        if td < 2.8:

            f = ((np.pi * td) ** -0.5 + 0.5 - (td / np.pi) ** 0.5 / 4 + td / 8)

        else:

            gamma = 0.57722
            theta = (np.log(4 * td) - 2 * gamma)
            f = 2 * (1 / theta - gamma / theta ** 2)

        return f

first_key = list(results.keys())[0]
time_array = np.array(results[first_key]["time"])
f_array = calculate_q(time_array)
f_array = np.array(f_array) * DT * 2 * np.pi

for key in results.keys():

    plt.plot(results[key]["time"], results[key]["value"], label=key)

plt.plot(results[key]["time"], f_array, "--", label="correlation", color="Black")

plt.xlabel("$t_d$ [-]")
plt.ylabel("$dT/dr$ [°C/m]")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
