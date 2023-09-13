from main_code.well_model.ground_models.FreeFemAnalysis.FreeFemAnalyzer import FreeFEMOptions
import os.path

from .edp_patches import (

    ProblemDefinitionOptions,
    MeshOptions,
    BASE_TXT

)


class TimeDependentOptions(FreeFEMOptions):

    def __init__(

        self, mesh_options: MeshOptions, problem_definition_options:ProblemDefinitionOptions, mesh_only=False

    ):

        self.mesh_options = mesh_options
        self.problem_definition_options = problem_definition_options
        self.mesh_only = mesh_only
        super().__init__()

    def edp_text(self):

        if not self.mesh_only:

            return BASE_TXT.format(

                mesh_definition = self.mesh_options.edp_text(),
                problem_definition = self.problem_definition_options.edp_text()

            )

        else:

            return BASE_TXT.format(

                mesh_definition=self.mesh_options.edp_text(),
                problem_definition=""

            )

    def other_calculation_folder_implementation(self):

        self.problem_definition_options.save_path = self.res_file_path

        if not self.mesh_options.retrieve_mesh:
            self.mesh_options.mesh_path = os.path.join(self.calculation_folder, self.mesh_options.mesh_path_name)

