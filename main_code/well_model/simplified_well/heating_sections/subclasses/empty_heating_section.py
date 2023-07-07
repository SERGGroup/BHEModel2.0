from main_code.well_model.simplified_well.heating_sections.abstract_class import AbstractHeatingSection
from scipy.optimize import Bounds


class EmptyHeatingSection(AbstractHeatingSection):

    def update_HS(self):

        """
            Simply equals the output and input conditions

        """

        self.output_point.set_variable("T", self.input_point.get_variable("T"))
        self.output_point.set_variable("rho", self.input_point.get_variable("rho"))

    def get_c_well(self):
        return 0.

    @property
    def initial_guess(self):
        return []

    @property
    def optimization_bounds(self):
        """
            Parameters to be optimized are:

                - dt_rocks
                - n_wells
                - d_well
                - tkn_annulus

        """
        return Bounds([], [])

    def set_optimization_param(self, optimization_param):
        pass

    def additional_setup_data(self, data_frame: dict):

        data_frame["Calculation Options"].update({

            "heating section": {"value": "EmptyHeatingSection", "unit": None}

        })

        return data_frame
