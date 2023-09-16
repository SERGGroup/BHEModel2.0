# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.ground_models.FreeFemAnalysis import (

    FreeFEMAnalyzer, TimeDependentOptions, MeshOptions, ProblemDefinitionOptions

)
from main_code.support.other.excel_exporter import write_excel_sheet
from main_code.constants import CALCULATION_FOLDER, os
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np


# %%------------   DEFINE WORKING FOLDERS                 -----------------------------------------------------------> #
calculation_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "2 - Calculations",
    "res", "calculation folder"

)

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "2 - Calculations",
    "res", "results"

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
a = 4
DT = 40.

n_mesh = 3
m_mesh = 20
time_steps = 10000
n_save = time_steps
results = dict()
Pe_list = [0, 0.1, 1, 10]

for Pe in Pe_list:

    key = "Pe={}".format(Pe)

    mo = MeshOptions(

        n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
        mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
        retrieve_mesh=False, add_visualization_mesh=True

    )

    pdo = ProblemDefinitionOptions(

        grad_rock=0., DT=DT,
        time_range=[1e-3, 1000], time_steps=time_steps,
        n_plot=0, n_save=n_save, pe=Pe, vx=1

    )

    tdo = TimeDependentOptions(mesh_options=mo, problem_definition_options=pdo, mesh_only=False)
    ffa = FreeFEMAnalyzer(options=tdo, calculation_folder=calculation_folder)

    print(key)
    gross_result = ffa.calculate()
    print(key)

    time_list = list()
    value_list = list()

    for line in gross_result:

        line_split = line.split(";")
        time_list.append(float(line_split[0].strip()))
        value_list.append(float(line_split[1].strip()) / (DT * 2 * np.pi))

    results.update({key:{

        "time": time_list,
        "value": value_list

    }})


# %%------------   SAVE RESULTS                           -----------------------------------------------------------> #
file_path = os.path.join(result_folder, "2 - effect of convection.xlsx")
keys_list = list(results.keys())

dataframe = {

    'Time': {"unit": ["-"], "values": [results[keys_list[0]]["time"]]},
    'value': {"unit": list(), "values": list()},

}
for Pe in Pe_list:

    key = "Pe={}".format(Pe)
    dataframe["value"]["unit"].append(key)
    dataframe["value"]["values"].append(results[key]["value"])

write_excel_sheet(excel_path=file_path, sheet_name="results", data_frame=dataframe, overwrite="hard")


# %%------------   REMOVE ERRORS                          -----------------------------------------------------------> #
nans = np.empty(len(results['Pe=10']['value']))
nans[:] = np.nan
results['Pe=10']['value'][875:] = nans[875:]


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
keys_list = list(results.keys())
fig, (ax1, ax3) = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(10, 4)

ax1.set_ylabel("$f$ [-]")
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("$t_d$ [-]")

key = "Pe={}".format(Pe_list[0])
time_array = results[key]["time"]
base_value_array = np.array(results[key]["value"])
ax1.plot(time_array, base_value_array, "--", color="black", alpha=0.6, linewidth=2, label=key)

steady_state_values = list()
for Pe in Pe_list[1:]:

    key = "Pe={}".format(Pe)
    value_array = np.array(results[key]["value"])
    steady_state_values.append(value_array[-1])

    ax1.plot(time_array, value_array, label=key)

steady_state_values[-1] = value_array[874]

ax1.legend(fontsize="8", loc=3)
ax1.set_title("Effect of Natural Convection")

l_reg = LinearRegression().fit(

    np.array([np.log10(Pe_list[1:])]).T, np.array([np.log10(steady_state_values)]).T

)

ax3.plot(Pe_list[1:], steady_state_values, "-x")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_ylabel("$f$ [-]")
ax3.set_xlabel("$Pe$ [-]")
ax3.set_title("Steady State Condition")

props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax3.text(

    0.1, 0.925,
    "$df\\ =\\ {C}\\ Pe^{{{a}}}$".format(

        a="{:.2f}".format(10 ** l_reg.coef_[0][0]),
        C="{:.2e}".format(10 ** l_reg.intercept_[0])

    ),
    verticalalignment='top', bbox=props,
    fontsize = 14, transform=ax3.transAxes

)

plt.tight_layout(pad=2)
plt.show()
