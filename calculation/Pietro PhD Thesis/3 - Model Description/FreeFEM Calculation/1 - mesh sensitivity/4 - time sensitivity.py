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

    CALCULATION_FOLDER, "Pietro PhD Thesis", "3 - Model Description",
    "FreeFEM Calculation", "1 - mesh sensitivity",
    "res", "calculation folder"

)

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis", "3 - Model Description",
    "FreeFEM Calculation", "1 - mesh sensitivity",
    "res", "results"

)


# %%------------   CALCULATIONS                           -----------------------------------------------------------> #
a = 4
DT = 40.
n_mesh = 2
m_mesh = 15
n_save = 500
results = dict()

for time_steps in [10, 50, 100]:

    print("{} time steps".format(time_steps))

    mo = MeshOptions(

        n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
        mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
        retrieve_mesh=False, add_visualization_mesh=True

    )

    if n_save > time_steps:
        n_save_real = time_steps

    else:
        n_save_real = n_save

    pdo = ProblemDefinitionOptions(

        grad_rock=0.01, DT=DT,
        time_range=[1e-3, 1000], time_steps=time_steps,
        n_plot=0, n_save=n_save_real

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

        "{} time steps".format(time_steps): {

            "time": time_list,
            "value": value_list

        }

    })

    print("{} time steps".format(time_steps))


# %%------------   SAVE RESULTS                           -----------------------------------------------------------> #
file_path = os.path.join(result_folder, "time-sensitivity.xlsx")

for key in results.keys():

    dataframe = {

        'Time': {"unit": ["-"], "values": [results[key]["time"]]},
        'value': {"unit": ["-"], "values": [results[key]["value"]]},

    }
    write_excel_sheet(excel_path=file_path, sheet_name=key, data_frame=dataframe, overwrite="hard")


# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
keys_list = list(results.keys())

keys_list.reverse()
for key in keys_list:

    plt.plot(results[key]["time"], results[key]["value"], label=key)

plt.xlabel("$t_d$ [-]")
plt.ylabel("$f$ [-]")
plt.xscale("log")
plt.yscale("log")
plt.legend(fontsize="8")
plt.title("Time Step Sensitivity")
plt.show()
