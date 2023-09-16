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
a = 4
n_mesh = 2
m_mesh = 20

mo = MeshOptions(

    n_points=a * n_mesh, n_points_circle=a * m_mesh * n_mesh, graph_r_ratio=15,
    mesh_path=os.path.join(calculation_folder, "mesh_out.mesh"),
    retrieve_mesh=False, add_visualization_mesh=True

)

pdo = ProblemDefinitionOptions(

    grad_rock=0.01, DT=40.,
    time_range=[0, 157.68], time_steps=10000

)

tdo = TimeDependentOptions(mesh_options=mo, problem_definition_options=pdo, mesh_only=True)
ffa = FreeFEMAnalyzer(options=tdo, calculation_folder=calculation_folder)
ffa.calculate()
