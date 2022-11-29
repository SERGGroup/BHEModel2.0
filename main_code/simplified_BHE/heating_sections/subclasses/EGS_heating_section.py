import numpy as np

from main_code.simplified_BHE.heating_sections.abstract_class import AbstractHeatingSection
from scipy.optimize import Bounds


class EGSHeatingSection(AbstractHeatingSection):

    def __init__(self, main_BHE, S_res=1.18e-12, k_res=5e-14):

        """

            Default values has been taken from:

            Adamns et Al.,
            A comparison of electric power output of CO2 Plume Geothermal (CPG) and brine
            geothermal systems for varying reservoir conditions, Applied Energy

        """

        super().__init__(main_BHE)

        self.S_res = S_res              # [m/s]     Average Specific Kinematic Viscosity
        self.k_res = k_res              # [m^2]     Horizontal Permeability
        # self.R_res = S_res / k_res      # [1/(m*s)] Average Specific Inverse Mobility

        self.R_res = 25

        self.dp_res = 0.

    def update_HS(self):

        self.dp_res = self.R_res * self.input_point.m_dot / 1e3    # Pa -> MPa
        self.output_point.set_variable("P", self.input_point.get_variable("P") - self.dp_res)
        self.output_point.set_variable("T", self.main_BHE.T_rocks)
        self.output_point.m_dot = self.input_point.m_dot

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

                - DT_rocks
                - n_wells
                - d_well
                - tkn_annulus

        """
        return Bounds([], [])

    def set_optimization_param(self, optimization_param):

        pass
