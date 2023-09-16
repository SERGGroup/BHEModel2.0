# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.ground_models.FreeFemAnalysis import (

    FreeFEMAnalyzer, TimeDependentOptions, MeshOptions, ProblemDefinitionOptions

)
from main_code.support.other.excel_exporter import write_excel_sheet
from main_code.constants import CALCULATION_FOLDER, os
from sklearn.linear_model import LinearRegression
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
DT = 40.
n_steps = 60
n_mesh_list = [1, 2, 3, 4]
m_mesh_list = [5, 10, 15, 20]

results = dict()

for m_mesh in m_mesh_list:

    for n_mesh in n_mesh_list:

        curr_key = "n={} - m={}".format(n_mesh, m_mesh)
        print(curr_key)

        mo = MeshOptions(

            n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
            mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
            retrieve_mesh=False, add_visualization_mesh=True

        )

        pdo = ProblemDefinitionOptions(

            grad_rock=0.01, DT=DT,
            time_range=[1e-3, 1000], time_steps=n_steps,
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
            value_list.append(float(line_split[1].strip()) / (DT * 2 * np.pi))

        results.update({

            curr_key: {

                "time": time_list,
                "value": value_list,
                "calc_time": ffa.mean_calculation_time

            }

        })

        print(curr_key)


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

    for j in range(1, len(res_keys_list)):

        if j > i:

            key_j = res_keys_list[j]
            arr_j = np.array(results[key_j]["value"])
            log_arr_j = np.log10(arr_j)

            # Standard Differences
            d_arr_sqr = (arr_j - arr_i) ** 2
            rel_diff = np.sqrt(np.mean(d_arr_sqr / (arr_i ** 2)))
            d_final = np.abs((arr_j[-1] - arr_i[-1])/arr_i[-1])


            # Log Based Differences
            log_d_arr_sqr = (log_arr_j - log_arr_i) ** 2
            log_rel_diff = np.sqrt(np.mean(log_d_arr_sqr / (log_arr_i ** 2)))
            log_d_final = np.abs((log_arr_j[-1] - log_arr_i[-1]) / log_arr_i[-1])

        else:

            d_final = "-"
            rel_diff = "-"
            log_d_final = "-"
            log_rel_diff = "-"

        difference_dataframe['relative difference']["values"][-1].append(rel_diff)
        difference_dataframe['final difference']["values"][-1].append(d_final)
        log_difference_dataframe['relative difference']["values"][-1].append(log_rel_diff)
        log_difference_dataframe['final difference']["values"][-1].append(log_rel_diff)


# %%------------   SAVE RESULTS                           -----------------------------------------------------------> #
file_path = os.path.join(result_folder, "mesh-sensitivity.xlsx")
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

for key in results.keys():

    plt.plot(results[key]["time"], results[key]["value"], label=key)

# plt.plot(time_array, f_array, "--", label="correlation", color="Black")

plt.xlabel("$t_d$ [-]")
plt.ylabel("$f$ [-]")
plt.xscale("log")
plt.yscale("log")
plt.legend(ncol=4, fontsize="8")
plt.title("Mesh Sensitivity")
plt.show()


# %%------------   PLOT ADDITIONAL RESULTS                -----------------------------------------------------------> #
last_key = "n={} - m={}".format(max(n_mesh_list), max(m_mesh_list))
fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(10, 4)

for m_mesh in m_mesh_list:

    comp_time_list = list()
    last_diff_list = list()
    diff_list = list()
    points_label=list()

    for n_mesh in n_mesh_list:

        curr_key = "n={} - m={}".format(n_mesh, m_mesh)

        if not curr_key == last_key:

            comp_time_list.append(results[curr_key]["calc_time"] / n_steps)

            arr_i = np.array(results[curr_key]["value"])
            arr_j = np.array(results[last_key]["value"])

            d_arr_sqr = (arr_j - arr_i) ** 2
            diff_list.append(np.sqrt(np.mean(d_arr_sqr / (arr_i ** 2))))
            last_diff_list.append(np.abs((arr_j[-1] - arr_i[-1]) / arr_i[-1]))

            points_label.append("n={}".format(n_mesh))

    ax1.plot(comp_time_list, diff_list, "-x", label="m={}".format(m_mesh))
    ax2.plot(comp_time_list, last_diff_list, "-x", label="m={}".format(m_mesh))

    if m_mesh == m_mesh_list[0]:

        for i in range(len(points_label)):

            ax1.text(comp_time_list[i]*1.1, diff_list[i], "${}$".format(points_label[i]), color="black", alpha=0.6)
            ax2.text(comp_time_list[i]*1.1, last_diff_list[i], "${}$".format(points_label[i]), color="black", alpha=0.6)

for ax in [ax1, ax2]:

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Time Step Evaluation Time [s/step]")
    ax.set_ylabel("Relative error [-]")
    ax.legend()

for n_mesh in n_mesh_list:

    comp_time_list = list()
    last_diff_list = list()
    diff_list = list()
    points_label=list()

    for m_mesh in m_mesh_list:

        curr_key = "n={} - m={}".format(n_mesh, m_mesh)

        if not curr_key == last_key:

            comp_time_list.append(results[curr_key]["calc_time"] / n_steps)

            arr_i = np.array(results[curr_key]["value"])
            arr_j = np.array(results[last_key]["value"])

            d_arr_sqr = (arr_j - arr_i) ** 2
            diff_list.append(np.sqrt(np.mean(d_arr_sqr / (arr_i ** 2))))
            last_diff_list.append(np.abs((arr_j[-1] - arr_i[-1]) / arr_i[-1]))

            points_label.append("n={}".format(n_mesh))

    ax1.plot(comp_time_list, diff_list, "--", color="black", alpha=0.6, linewidth=0.75)
    ax2.plot(comp_time_list, last_diff_list, "--", color="black", alpha=0.6, linewidth=0.75)

ax2.set_title("Last Time-Step Error")
ax1.set_title("Overall Error")
plt.tight_layout(pad=2)
plt.show()

# %%------------   ACCURACY LINEAR REGRESSION             -----------------------------------------------------------> #
x = np.zeros((len(m_mesh_list)*len(n_mesh_list) - 1, 2))
y = np.zeros((len(m_mesh_list)*len(n_mesh_list) - 1, 3))
last_key = "n={} - m={}".format(max(n_mesh_list), max(m_mesh_list))

i = 0
for m_mesh in m_mesh_list:

    for n_mesh in n_mesh_list:

        curr_key = "n={} - m={}".format(n_mesh, m_mesh)

        if not curr_key == last_key:

            x[i, 0] = m_mesh
            x[i, 1] = n_mesh
            y[i, 0] = results[curr_key]["calc_time"] / n_steps

            arr_i = np.array(results[curr_key]["value"])
            arr_j = np.array(results[last_key]["value"])
            d_arr_sqr = (arr_j - arr_i) ** 2
            y[i, 1] = np.sqrt(np.mean(d_arr_sqr / (arr_i ** 2)))
            y[i, 2] = np.abs((arr_j[-1] - arr_i[-1]) / arr_i[-1])

            i += 1

x = np.log10(x)
y = np.log10(y)
models = list()

for k in range(3):

    models.append(LinearRegression().fit(x, y[:,k]))
