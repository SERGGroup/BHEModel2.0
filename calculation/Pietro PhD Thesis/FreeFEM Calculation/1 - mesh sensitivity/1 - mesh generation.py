# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_code.well_model.ground_models.FreeFemAnalysis import (

    FreeFEMAnalyzer, TimeDependentOptions, MeshOptions, ProblemDefinitionOptions

)
from main_code.constants import CALCULATION_FOLDER, os


# %%------------   DEFINE WORKING FOLDERS                 -----------------------------------------------------------> #
calculation_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "FreeFEM Calculation", "1 - mesh sensitivity",
    "res", "calculation folder"

)


# %%------------   INITIALIZATION                         -----------------------------------------------------------> #
mo = MeshOptions(

    n_points=10, n_points_circle=30, graph_r_ratio=25,
    ratio_L=1500., ratio_H=500., ref_index=2,
    add_visualization_mesh=True, show_mesh=True

)

pdo = ProblemDefinitionOptions(

    grad_rock=0.01, DT=40.,
    max_time=157.68, time_steps=10000

)

tdo = TimeDependentOptions(mesh_options=mo, problem_definition_options=pdo, mesh_only=True)
ffa = FreeFEMAnalyzer(options=tdo, calculation_folder=calculation_folder)
ffa.calculate()
