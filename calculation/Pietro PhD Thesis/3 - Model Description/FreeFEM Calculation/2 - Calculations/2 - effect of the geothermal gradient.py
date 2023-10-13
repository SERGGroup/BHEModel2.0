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

    CALCULATION_FOLDER, "Pietro PhD Thesis", "3 - Model Description",
    "FreeFEM Calculation", "2 - Calculations",
    "res", "calculation folder"

)

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis", "3 - Model Description",
    "FreeFEM Calculation", "2 - Calculations",
    "res", "results"

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
a = 4
DT = 40.

n_mesh = 3
m_mesh = 20
time_steps = 1000
n_save = time_steps
results = dict()
nd_grads = [0, 0.1, 1, 10]

for nd_grad in nd_grads:

    key = "nd_grad={}".format(nd_grad)

    mo = MeshOptions(

        n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
        mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
        retrieve_mesh=False, add_visualization_mesh=True

    )

    pdo = ProblemDefinitionOptions(

        grad_rock=nd_grad*DT, DT=DT,
        time_range=[1e-3, 1000], time_steps=time_steps,
        n_plot=0, n_save=n_save

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
file_path = os.path.join(result_folder, "2 - gradient effect.xlsx")
keys_list = list(results.keys())

dataframe = {

    'Time': {"unit": ["-"], "values": [results[keys_list[0]]["time"]]},
    'value': {"unit": list(), "values": list()},

}
for nd_grad in nd_grads:

    key = "nd_grad={}".format(nd_grad)
    dataframe["value"]["unit"].append(key)
    dataframe["value"]["values"].append(results[key]["value"])

write_excel_sheet(excel_path=file_path, sheet_name="results", data_frame=dataframe, overwrite="hard")


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
keys_list = list(results.keys())
fig, (ax1, ax3) = plt.subplots(1, 2, dpi=300)
fig.set_size_inches(10, 4)

ax1.set_ylabel("$f$ [-]")
ax3.set_ylabel("$df$ [-]")

for ax in [ax1, ax3]:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$t_d$ [-]")

ax2 = ax1.inset_axes([0.5, 0.5, 0.4, 0.4])
x1, x2, y1, y2 = 1e2, 1e3, 2.45e-1, 3.5e-1
ax2.set_xlim(x1, x2)
ax2.set_ylim(y1, y2)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax1.indicate_inset_zoom(ax2, edgecolor="black")

key = "nd_grad={}".format(nd_grads[0])
time_array = results[key]["time"]
base_value_array = np.array(results[key]["value"])
ax1.plot(time_array, base_value_array, "--", color="black", alpha=0.6, linewidth=2, label=key)
ax2.plot(time_array, base_value_array, "--", color="black", alpha=0.6, linewidth=2)

last_values_diff = list()
for nd_grad in nd_grads[1:]:

    key = "nd_grad={}".format(nd_grad)
    value_array = np.array(results[key]["value"])
    diff_array = value_array - base_value_array
    last_values_diff.append(diff_array[-1])

    ax3.plot(time_array, diff_array, label=key)
    ax1.plot(time_array, value_array, label=key)
    ax2.plot(time_array, value_array)

ax1.legend(fontsize="8", loc=3)
ax3.legend(fontsize="8")
ax1.set_title("Effect of Geothermal Gradient")
ax3.set_title("Difference with No Gradient")
plt.tight_layout(pad=2)
plt.show()

# %%------------   ADDITIONAL PLOT                        -----------------------------------------------------------> #

l_reg = LinearRegression().fit(

    np.array([np.log10(nd_grads[1:])]).T, np.array([np.log10(last_values_diff)]).T

)

fig, ax = plt.subplots()

ax.plot(nd_grads[1:], last_values_diff, "-x")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("$df$ [-]")
ax.set_xlabel("${\\nabla T_{rocks}}_{d}$ [-]")
ax.set_title("Difference-Gradient Correlation")

props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(

    0.05, 0.925,
    "$df\\ =\\ {C}\\ {{{{\\nabla T_{{rocks}}}}_{{d}}}}$".format(

        C="{:.3e}".format(10 ** l_reg.intercept_[0])

    ),
    verticalalignment='top', bbox=props,
    fontsize = 14, transform=ax.transAxes
)

plt.show()
