from REFPROPConnector.Tools.units_converter import convert_variable
from .abstract_surface_plant import AbstractSurfacePlant
import numpy as np


class BaseSurfacePlant(AbstractSurfacePlant):

    def __init__(self, bhe_well):

        super().__init__()
        self.bhe_well = bhe_well

        self.w_dot = 0.
        self.ex_dot = 0.
        self.eta_exs = 0.

        self.w_dot_nd = np.nan
        self.ex_dot_nd = np.nan
        self.w_dex_min = np.nan
        self.w_dex_max = np.nan

    def prepare_calculation(self):

        self.bhe_well.update_simplified()
        self.append_BHE_well_points(

            BHE_input=self.bhe_well.points[0],
            BHE_output=self.bhe_well.points[-1]

        )

        self.append_ambient_condition(self.bhe_well.points[0].duplicate())

    def thermo_analysis(self):

        self.points.append(self.BHE_well_points[0])
        self.points.append(self.BHE_well_points[0].duplicate())
        self.points.append(self.BHE_well_points[0].duplicate())
        self.points.append(self.BHE_well_points[1])

        self.points[1].set_variable("P", self.points[0].get_variable("P"))
        self.points[1].set_variable("S", self.points[-1].get_variable("S"))

        self.points[2].set_variable("P", self.points[-1].get_variable("P"))
        self.points[2].set_variable("S", self.points[0].get_variable("S"))

        t_amb_k, info = convert_variable(self.points[-1].get_variable("T"), "T", self.points[-1].get_unit("T"), "K")

        self.w_dot = self.points[0].get_variable("H") - self.points[-1].get_variable("H")
        self.ex_dot = self.w_dot - t_amb_k * (self.points[0].get_variable("S") - self.points[-1].get_variable("S"))

        if (not (self.w_dot < 0 or self.ex_dot < 0)) and (not self.points[0].get_variable("H") < -1e6):

            cf = 1 -  t_amb_k / (self.bhe_well.t_rocks + 273.15)

            self.w_dot_nd = self.w_dot / (self.points[-1].get_variable("CP") * self.points[-1].get_variable("T"))
            self.ex_dot_nd = self.ex_dot / (self.points[-1].get_variable("CP") * self.points[-1].get_variable("T"))
            self.eta_exs = self.ex_dot / (self.w_dot * cf)

            self.w_dex_min = (self.points[1].get_variable("H") - self.points[-1].get_variable("H")) / self.w_dot
            self.w_dex_max = (self.points[0].get_variable("H") - self.points[2].get_variable("H")) / self.w_dot

    def economic_analysis(self):
        pass

    def LCA_analysis(self):
        pass