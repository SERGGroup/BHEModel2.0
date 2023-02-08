import os

CURRENT_FOLDER = os.path.dirname(os.path.abspath(__file__))
MAIN_FOLDER = os.path.abspath(os.path.join(CURRENT_FOLDER, os.pardir))
CALCULATION_FOLDER = os.path.join(MAIN_FOLDER, "calculation")
RES_FOLDER = os.path.join(CALCULATION_FOLDER, "resources")
TEST_RES_FOLDER = os.path.join(RES_FOLDER, "test")
EES_RES_FOLDER = os.path.join(RES_FOLDER, "EESfiles")
GROUND_MESH_FOLDER = os.path.join(RES_FOLDER, "GroundMesh")
FORTRAN_CODE_FOLDER = os.path.join(RES_FOLDER, "fortran_code")
FREE_FEM_FOLDER = os.path.join(RES_FOLDER, "FreeFEM_files")
OTHER_FOLDER = os.path.join(RES_FOLDER, "support")
