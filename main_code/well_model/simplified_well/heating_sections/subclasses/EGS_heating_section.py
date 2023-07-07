from main_code.well_model.simplified_well.heating_sections.abstract_class import AbstractHeatingSection
from scipy.optimize import Bounds
import numpy as np


class EGSHeatingSection(AbstractHeatingSection):

    def __init__(self, main_BHE, S_res=1.18e-12, k_res=5e-14, R_res=None, Q_res=None):

        """

            Default values has been taken from:

            Adamns et Al.,
            A comparison of electric power output of CO2 Plume Geothermal (CPG) and brine
            geothermal systems for varying reservoir conditions, Applied Energy

        """

        super().__init__(main_BHE)

        self.S_res = S_res              # [m/s]     Average Specific Kinematic Viscosity
        self.k_res = k_res              # [m^2]     Horizontal Permeability
        self.R_res = S_res / k_res      # [1/(m*s)] Average Specific Inverse Mobility

        if R_res is not None:

            self.R_res = R_res

        self.Q_res = Q_res
        self.dp_res = 0.

    def update_HS(self):

        self.dp_res = self.R_res * self.input_point.m_dot / 1e3    # Pa -> MPa
        self.output_point.m_dot = self.input_point.m_dot

        if self.Q_res is None:

            self.output_point.set_variable("P", self.input_point.get_variable("P") - self.dp_res)
            self.output_point.set_variable("T", self.main_BHE.t_rocks)

        else:

            h_out = self.input_point.get_variable("h") + self.Q_res / self.input_point.m_dot

            self.output_point.set_variable("P", self.input_point.get_variable("P") - self.dp_res)
            self.output_point.set_variable("h", h_out)

    def get_c_well(self):

        #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #   !!!!!  AS NO HEATING SECTION GEOMETRY HAS BEEN DEFINED,   !!!!!
        #   !!!!!  THE CALCULATION ARE BASED ON THE WELL DEPTH ONLY   !!!!!
        #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        d_well = 0.254
        self.l_HS = 15 * self.main_BHE.dh * self.main_BHE.m_dot
        l_tot = self.main_BHE.dz_well + self.l_HS

        c_well = (0.105 * l_tot ** 2 + 1776 * l_tot * d_well + 2.753 * pow(10, 5))

        if self.main_BHE.working_fluid == "CarbonDioxide":

            dc_CO2 = (133 + 256 * d_well) * l_tot

        else:

            dc_CO2 = 0.

        return c_well, dc_CO2

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

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

            "heating section": {"value": "EGSHeatingSection", "unit": None}

        })

        data_frame.update({"Heating Section Data": {

            "S_res": {"value": self.S_res, "unit": "[m/s]"},
            "k_res": {"value": self.k_res, "unit": "[m^2]"},
            "R_res": {"value": self.R_res, "unit": "[1/(m*s)]"},

        }})

        return data_frame
