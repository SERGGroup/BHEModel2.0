from REFPROPConnector.Tools.units_converter import convert_variable
from .abstract_surface_plant import AbstractSurfacePlant
from ..support import PlantThermoPoint
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

        if (not (self.w_dot <= 0 or self.ex_dot <= 0)) and (not self.points[0].get_variable("H") < -1e6):

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


class BaseCO2HTHP(AbstractSurfacePlant):

    def __init__(self, bhe_well, p_steam, q_rel=0, p_max=50):

        super().__init__()
        self.bhe_well = bhe_well
        self.p_steam = p_steam
        self.q_rel = q_rel
        self.p_max = p_max

        self.m_steam = 1
        self.eta_comp = 0.80
        self.eta_turb = 0.75
        self.eta_pump = 0.8
        self.dt_he = 8

        self.is_p_max_limited = False
        self.is_physical = True

    def __init_points(self):

        self.points = list()

        # Point 0 - BHE Output
        self.points.append(self.BHE_output)

        # Points 1-3 - CO2 BHE 2022-10-04 - HTHP internal points:
        #
        #   1 - CO2 Compressor Output
        #   2 - CO2 Main HE Output
        #   3 - CO2 Expander Output

        for i in range(3):
            self.points.append(PlantThermoPoint(["CarbonDioxide"], [1]))

        # Point 4 - BHE Input

        self.points.append(self.BHE_input)

        # Points 5-9 - Steam Drum internal points:
        #
        #   5 - Steam Drum Liquid Output
        #   6 - Circulation Pump Output
        #   7 - Water Main HE Output
        #   8 - Lamination Valve Output
        #   9 - Steam Drum Vapour Output

        for i in range(5):
            self.points.append(PlantThermoPoint(["water"], [1]))

        # point 5 (saturated liquid) and 9 (saturated vapour) properties will remain fixed (only their flow rates
        # will change)

        self.points[5].set_variable("P", self.p_steam)
        self.points[5].set_variable("x", 0)

        self.points[9].set_variable("P", self.p_steam)
        self.points[9].set_variable("x", 1)

    def prepare_calculation(self):

        self.bhe_well.update_simplified()
        self.append_BHE_well_points(

            BHE_input=self.bhe_well.points[0],
            BHE_output=self.bhe_well.points[-1]

        )

        self.append_ambient_condition(self.bhe_well.points[0].duplicate())
        self.__init_points()

    def thermo_analysis(self):

        # STEP - 1
        # Calculate point 3 given q_rel
        self.points[3].set_variable("P", self.points[4].get_variable("P"))
        self.points[3].set_variable("H", self.points[4].get_variable("H") - self.q_rel)

        # Point 6 - Circulation Pump Output
        p_water = self.points[5].RPHandler.PC
        self.points[6].set_to_compression_result(p_water, self.eta_pump, input_state=self.points[5])

        # STEP - 2
        # Iterate point 2 and retrieve pressure
        p_up, t_up = self.__iterate_p_up()
        self.points[2].set_variable("P", p_up)
        self.points[2].set_variable("T", t_up)

        # STEP - 3
        # Evaluate point 1
        self.points[1].set_to_compression_result(p_up, self.eta_comp, self.points[0])

        # Point 7 - Water Main HE Output
        self.points[7].set_variable("T", self.points[1].get_variable("T") - self.dt_he)
        self.points[7].set_variable("P", p_water)

        # Point 8 - Lamination Valve Output
        self.points[8].set_variable("h", self.points[7].get_variable("h"))
        self.points[8].set_variable("P", self.p_steam)

        self.is_physical = self.points[6].get_variable("T") < self.points[2].get_variable("T")
        self.__evaluate_flow_rates()
        self.__calculate_power()

    def __evaluate_flow_rates(self):

        if 0 < self.points[8].get_variable("x") < 1:
            steam_water_flow_ratio = self.points[8].get_variable("x")

        elif self.points[8].get_variable("T") > self.points[5].get_variable("T"):
            steam_water_flow_ratio = 1

        else:
            steam_water_flow_ratio = -1

        dh_CO2 = self.points[1].get_variable("h") - self.points[2].get_variable("h")
        dh_water = self.points[7].get_variable("h") - self.points[6].get_variable("h")
        self.steam_CO2_flow_ratio = dh_CO2 / dh_water * steam_water_flow_ratio

        if self.steam_CO2_flow_ratio == 0:

            self.m_BHE = 99999999999999999

        else:

            self.m_BHE = self.m_steam / self.steam_CO2_flow_ratio

        self.m_dot_ratio_real = self.m_steam / self.m_BHE
        self.m_dot_ratio = 2.45 * self.m_steam / self.m_BHE

        # set flow rates to points

        # CO2 points
        for i in range(5):
            self.points[i].m_dot = self.m_BHE

        # Water points
        for i in range(5, 9):
            self.points[i].m_dot = self.m_steam / steam_water_flow_ratio

        # Steam point
        self.points[9].m_dot = self.m_steam

    def __calculate_power(self):

        m_dot_co2 = self.points[0].get_variable("m_dot")
        m_dot_water = self.points[5].get_variable("m_dot")

        # High Temperature Heat Exchanger
        dh_HTHE = self.points[1].get_variable("h") - self.points[2].get_variable("h")
        self.Q_HTHE = dh_HTHE * m_dot_co2

        # Gas Cooler
        dh_GC = self.points[3].get_variable("h") - self.points[4].get_variable("h")
        self.Q_GC = dh_GC * m_dot_co2

        # Turbine
        dh_turb = self.points[2].get_variable("h") - self.points[3].get_variable("h")
        self.W_turb = m_dot_co2 * dh_turb

        # Compressor
        dh_comp = self.points[1].get_variable("h") - self.points[0].get_variable("h")
        self.W_comp = m_dot_co2 * dh_comp

        # pump
        dh_pump = self.points[6].get_variable("h") - self.points[5].get_variable("h")
        self.W_pump = m_dot_water * dh_pump

        # BHE
        dh_BHE = self.BHE_output.get_variable("h") - self.BHE_input.get_variable("h")
        self.Q_BHE = m_dot_co2 * dh_BHE

        self.Q_ratio = self.Q_HTHE / self.Q_BHE

        if (self.W_comp + self.W_pump) > self.W_turb:

            self.COP = self.Q_HTHE / (self.W_comp + self.W_pump - self.W_turb)

        else:

            self.COP = (self.Q_HTHE + self.W_turb - (self.W_comp + self.W_pump)) / (self.W_comp + self.W_pump)

    def __iterate_p_up(self):

        # iteration uses bisection method for speed, P_up initial is P_CO2_max P_down initial is P_BHE_out
        t_2 = self.points[6].get_variable("T") + self.dt_he
        t_goal = self.points[3].get_variable("T")

        if t_2 < t_goal:

            self.is_p_max_limited = False
            return self.points[3].get_variable("P")

        toll = 0.001
        error = 10 * toll
        self.iter = 0

        p_out = self.points[3].get_variable("P")

        tmp_point_in = self.points[2].duplicate()
        tmp_point_out = self.points[2].duplicate()

        tmp_point_in.set_variable("T", t_2)
        tmp_point_in.set_variable("P", self.p_max)
        tmp_point_out.set_to_expansion_result(p_out, self.eta_turb, tmp_point_in)

        if tmp_point_out.get_variable("T") > t_goal:

            self.is_p_max_limited = True
            p_limit = [self.p_max, self.p_max]
            t_limit = [t_goal, t_2]

        else:

            self.is_p_max_limited = False
            p_limit = [p_out, self.p_max]
            t_limit = [t_2, t_2]

        p_mean = sum(p_limit) / 2
        t_mean = sum(t_limit) / 2

        while error > toll:

            tmp_point_in.set_variable("P", p_mean)
            tmp_point_in.set_variable("T", t_mean)
            tmp_point_out.set_to_expansion_result(p_out, self.eta_turb, tmp_point_in)

            t_res = tmp_point_out.get_variable("T")

            if t_res < t_goal:

                t_limit[0] = t_mean
                p_limit[1] = p_mean

            else:

                t_limit[1] = t_mean
                p_limit[0] = p_mean

            error = abs(t_res - t_goal)
            p_mean = sum(p_limit) / 2
            t_mean = sum(t_limit) / 2

            self.iter += 1

            if self.iter > 200:
                break

        return p_mean, t_mean

    def economic_analysis(self):
        pass

    def LCA_analysis(self):
        pass