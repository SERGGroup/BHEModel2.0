from main_code.simplified_well.heating_sections.abstract_class import AbstractHeatingSection
from scipy.optimize import Bounds


class DefaultHeatingSection(AbstractHeatingSection):

    def update_HS(self):
        self.output_point.set_variable("P", self.input_point.get_variable("P"))
        self.output_point.set_variable("T", self.main_BHE.T_rocks)

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