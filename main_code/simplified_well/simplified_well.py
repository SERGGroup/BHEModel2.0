from main_code.simplified_well.heating_sections import DefaultHeatingSection, AbstractHeatingSection
from EEETools.Tools.API.ExcelAPI.modules_importer import calculate_excel
from main_code.support.other.simple_integrator import SimpleIntegrator
from EEETools.Tools.API.Tools.main_tools import get_result_data_frames
from main_code.support import PlantThermoPoint, retrieve_PPI
from abc import ABC, abstractmethod
from scipy.integrate import RK45
from scipy.optimize import Bounds
from main_code import constants
from scipy.constants import g
import numpy as np
from sty import ef
import warnings
import os


class SimplifiedWell(ABC):

    def __init__(

            self, input_thermo_point: PlantThermoPoint, dz_well, T_rocks,
            heating_section=None, k_rocks=0.2, c_rocks=1, rho_rocks=2500,
            t_surf=10, geo_flux=0.1, q_up_down=0.0, PPI=None, use_rk=True,
            d_inj=None, d_prod=None, discretize_p_losses=False

    ):

        # GEOTHERMAL FIELD CONDITION

        self.dz_well = dz_well      # unit: m .............. default: NONE (REQUIRED)
        self.T_rocks = T_rocks      # unit: °C ............. default: NONE (REQUIRED)
        self.geo_flux = geo_flux    # unit: W / m^2 ........ default: 0.1 W / m^2

        self.k_rocks = k_rocks      # unit: W / (m K) ...... default: 0.2 W / (m K)
        self.c_rocks = c_rocks      # unit: kJ / (kg K) .... default: 1 kJ / (kg K)
        self.rho_rocks = rho_rocks  # unit: kg / m^3 ....... default: 2500 kg / m^3

        self.d_inj = d_inj          # unit: m .............. default: NONE (PRESSURE LOSSES IGNORED)
        self.d_prod = d_prod        # unit: m .............. default: NONE (PRESSURE LOSSES IGNORED)

        #   alpha_rocks -> rock thermal diffusivity in [m^2/s]
        #   (1e3 conversion factor c_rocks [kJ / (kg K)] -> [J / (kg K)])
        self.alpha_rocks = k_rocks / (rho_rocks * c_rocks * 1e3)
        self.ovr_grad = (T_rocks - t_surf) / dz_well

        self.geo_gradient = None
        self.depth_optimization = False
        self.use_rk = use_rk
        self.discretize_p_losses = discretize_p_losses

        self.q_dot_up = - q_up_down
        self.q_dot_down = q_up_down

        self.c_well = 0.
        self.C0 = [0., 0.]
        self.P_loss = [0., 0.]

        self.integrators_profiler = list()
        self.__init_points(input_thermo_point)
        self.__init_heating_section(heating_section)
        self.__reset_control_elements(first_initialization=True)
        self.__ambient_condition = None

        self.__current_PPI = PPI

    def __init_points(self, input_thermo_point):

        self.points = [input_thermo_point]
        self.__tmp_point = input_thermo_point.duplicate()

    def __init_heating_section(self, heating_section):

        self.__heating_section_changed = False
        self.heating_section = heating_section

    def __reset_control_elements(self, first_initialization=False):

        self.__old_input_values = {

            "T": self.points[0].get_variable("T"),
            "P": self.points[0].get_variable("P"),
            "m_dot": self.points[0].get_variable("m_dot")

        }

        if first_initialization:
            self.__heating_section_changed = True

        else:
            self.__heating_section_changed = False

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------   PROPERTIES   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def current_PPI(self):

        if self.__current_PPI is None:

            self.__current_PPI = retrieve_PPI()

        return self.__current_PPI

    @current_PPI.setter
    def current_PPI(self, PPI):

        self.__current_PPI = PPI

    @property
    def m_dot(self):

        return self.points[0].m_dot

    @property
    def l_tot(self):

        return self.dz_well + self.heating_section.l_HS

    @property
    def conditions_changed(self):

        try:

            input_cond_changed = True
            for param in ["T", "P", "m_dot"]:
                input_cond_changed = input_cond_changed and self.__old_input_values[param] == self.points[0].get_variable(param)

            if input_cond_changed and self.__heating_section_changed:
                return True

        except:

            return True

        return self.heating_section.condition_changed

    @property
    def working_fluid(self):

        try:

            working_fluid = self.points[0].RPHandler.points[0]

        except:

            working_fluid = ""

        return working_fluid

    @property
    def heating_section(self):

        return self.__heating_section

    @heating_section.setter
    def heating_section(self, input_heating_section):

        if not issubclass(type(input_heating_section), AbstractHeatingSection):

            input_heating_section = DefaultHeatingSection(self)

        self.__heating_section = input_heating_section

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update(self, ambient_condition=None):

        self.__set_reference_state(ambient_condition)

        self.integrators_profiler = list()
        self.evaluate_points()

        self.__evaluate_parameters()
        self.__evaluate_performance_parameters()
        self.__perform_exergy_analysis()

        self.__reset_control_elements()

    def update_simplified(self):

        self.integrators_profiler = list()
        self.evaluate_points()

    def __set_reference_state(self, ambient_condition):

        self.__ambient_condition = ambient_condition

        if type(ambient_condition) == PlantThermoPoint:

            ref_state = ambient_condition

        else:

            ref_state = self.points[0].duplicate()

        self.points[0].reference_state = ref_state

    @abstractmethod
    def evaluate_points(self):

        """

            METHOD TO BE IMPLEMENTED IN SUBCLASSES
            evaluate the well points according to the reservoir modelling hypothesis

        """

    def __evaluate_parameters(self):

        P_list = list()
        T_list = list()
        h_list = list()
        s_list = list()

        for point in self.points:

            P_list.append(point.get_variable("P"))
            T_list.append(point.get_variable("T"))
            h_list.append(point.get_variable("h"))
            s_list.append(point.get_variable("s"))

        self.P_ratio = P_list[3] / P_list[0]

        self.dP = P_list[3] - P_list[0]
        self.dT = T_list[3] - T_list[0]
        self.dh = h_list[3] - h_list[0]
        self.ds = s_list[3] - s_list[0]

        self.dex = self.dh - self.points[0].RPHandler.T_0_in_K * self.ds
        self.q_bottom = h_list[2] - h_list[1]
        self.Q_bottom = self.q_bottom * self.points[0].m_dot
        self.power = self.dh * self.points[0].m_dot

        self.eta_I = self.dh / self.q_bottom

    def __evaluate_performance_parameters(self):

        P_list = list()
        h_list = list()
        s_list = list()
        rho_list = list()

        for point in self.points:
            P_list.append(point.get_variable("P"))
            h_list.append(point.get_variable("h"))
            s_list.append(point.get_variable("s"))
            rho_list.append(point.get_variable("rho"))

        max_dP = rho_list[0] * g * self.dz_well / 10 ** 6
        self.eta_pressure = self.dP / max_dP
        self.beta_turb = P_list[1] / P_list[0]
        self.P_ratio_depth = self.P_ratio / self.dz_well
        self.DP_depth = self.dP / self.dz_well
        self.Dh_depth = self.dh / self.dz_well
        self.DT_depth = self.dT / self.dz_well

    def __perform_exergy_analysis(self):

        excel_path = os.path.join(constants.RES_FOLDER, "3ETool_Excel-Files", "BHE_simple_analysis.xlsx")
        T_rocks_K = self.T_rocks + 273.15

        new_exergy_list = [

            {"index": 20, "value": self.dz_well * g / 1000},
            {"index": 21, "value": (1 - self.points[0].RPHandler.T_0_in_K / T_rocks_K) * self.q_bottom},
            {"index": 22, "value": self.dex}

        ]

        for index in range(len(self.points)):
            new_exergy_list.append({

                "index": index + 1,
                "value": self.points[index].get_variable("exergy")

            })

        with warnings.catch_warnings():

            warnings.simplefilter("ignore")

            array_handler = calculate_excel(excel_path, new_exergy_list=new_exergy_list, export_solution=False)
            result = get_result_data_frames(array_handler)

        self.comp_results = result["Comp Out"]
        self.eta_II = array_handler.overall_efficiency

    def update_DP_vertical(self, input_point: PlantThermoPoint, is_upward=True):

        if input_point == self.points[-1]:

            self.__tmp_point = input_point.duplicate()
            self.points.append(self.__tmp_point)

        else:

            self.__tmp_point = self.points[self.points.index(input_point) + 1]

        if is_upward:

            self.__mult = -1
            dq = self.q_dot_up

        else:

            self.__mult = 1
            dq = self.q_dot_down

        rho_in = input_point.get_variable("rho")
        dh = (self.__mult * g / 1e3 + dq) * self.dz_well
        self.__dh_dp_stream = (g + dq * 1e3) / g * 1e3

        p0 = input_point.get_variable("P")
        rho0 = input_point.get_variable("rho")
        self.integrators_profiler.append(list())

        if not self.use_rk:

            integrator = SimpleIntegrator(self.rk_funct, 0, [p0, rho0], self.dz_well, n_steps=100)

        else:

            integrator = RK45(self.rk_funct, 0, [p0, rho0], self.dz_well)

        while integrator.status == 'running':

            integrator.step()

            range_list = [integrator.t_old, integrator.t]
            range_list.sort()

            self.integrators_profiler[-1].append({

                "range": range_list,
                "dense_out": integrator.dense_output(),
                "error": (not integrator.status == 'failed')

            })

        output = integrator.y
        self.__tmp_point.set_variable("h", input_point.get_variable("h") + dh)
        self.__tmp_point.set_variable("P", output[0])

        co_in = self.__calculate_C0(input_point, self.__dh_dp_stream)
        co_out = self.__calculate_C0(self.__tmp_point, self.__dh_dp_stream)
        c0 = (co_in + co_out) / 2

        if not self.discretize_p_losses:

            dp_loss = self.__evaluate_pressure_losses(input_point, self.__tmp_point, c0)
            self.__tmp_point.set_variable("P", self.__tmp_point.get_variable("P") - dp_loss)
            self.__tmp_point.set_variable("h", input_point.get_variable("h") + dh)

        else:

            dp_loss = 0.

        self.__tmp_point = self.__tmp_point.duplicate()

        return c0, dp_loss

    def rk_funct(self, z, y):

        p_curr = y[0]
        rho_curr = y[1]

        # h_curr = input_point.get_variable("h") + (g / 1e3 + dq) * z

        self.__tmp_point.set_variable("P", p_curr)
        self.__tmp_point.set_variable("rho", rho_curr)

        c0_curr = self.__calculate_C0(self.__tmp_point, self.__dh_dp_stream)

        if self.discretize_p_losses:

            dp_loss = self.__evaluate_pressure_losses(self.points[0], self.__tmp_point, dz=1)

        else:

            dp_loss = 0.

        dp = self.__mult * rho_curr * 9.81 / 1e6 - dp_loss
        d_rho = c0_curr * dp

        return [dp, d_rho]

    def __evaluate_pressure_losses(self, input_point, output_point, c0=None, dz=None):

        m_dot = input_point.m_dot
        g = 9.81

        if dz is None:
            dz = abs(self.dz_well)

        if input_point.get_variable("P") > output_point.get_variable("P"):
            d = self.d_prod

        else:
            d = self.d_inj

        if d is not None:

            mu = (input_point.get_variable("mu") + output_point.get_variable("mu")) / 2
            re = 4 * m_dot / (np.pi * d * mu / 1e6)

            f = self.__calculate_friction_factor(d, re)

            if c0 is not None:

                rho0 = max(input_point.get_variable("rho"), output_point.get_variable("rho"))
                l_mod = (np.exp(c0 / 1e6 * g * dz) - 1) / (g * c0 / 1e6)
                dp = 8 * (f * m_dot ** 2) / (d ** 5 * np.pi ** 2 * rho0) * l_mod

            else:

                rho0 = output_point.get_variable("rho")
                dp = 8 * (f * m_dot ** 2) / (d ** 5 * np.pi ** 2 * rho0) * dz

        else:

            dp = 0.

        return dp / 1e6

    @staticmethod
    def __calculate_friction_factor(d, re):

        if d is not None:

            rough = 55 * 1e-6       # [m] surface roughness
            rel_rough = rough / d

            A = (2.457 * np.log(1. / ((7. / re) ** 0.9 + 0.27 * rel_rough))) ** 16
            B = (37530. / re) ** 16

            f = 8. * ((8. / re) ** 12 + (1 / (A + B)) ** 1.5) ** (1. / 12.)

        else:

            f = 0.

        return f

    @staticmethod
    def __calculate_C0(thermo_point: PlantThermoPoint, dh_dp_stream):

        d_rho_d_P = thermo_point.get_derivative("rho", "P", "T")
        d_rho_d_T = thermo_point.get_derivative("rho", "T", "P")
        d_h_d_T = thermo_point.get_derivative("h", "T", "P")
        d_h_d_P = thermo_point.get_derivative("h", "P", "T")
        rho = thermo_point.get_variable("rho")

        return d_rho_d_P + d_rho_d_T / d_h_d_T * (dh_dp_stream / rho - d_h_d_P)

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update_economic(self, s_well=0.95):

        self.s_well = s_well
        if self.conditions_changed:

            self.update(self.__ambient_condition)

        self.heating_section.evaluate_cost()
        self.c_well = self.heating_section.CAPEX_well

    def force_PPI_update(self):

        self.current_PPI = None
        return self.current_PPI

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def initial_guess(self):

        if self.geo_gradient is not None:

            initial_guess = [1000]
            self.depth_optimization = True

        else:

            initial_guess = []
            self.depth_optimization = False

        initial_guess.extend(self.heating_section.initial_guess)
        return initial_guess

    @property
    def optimization_bounds(self):

        if self.geo_gradient is not None:

            lb = [500]
            ub = [5000]
            self.depth_optimization = True

        else:

            lb = []
            ub = []
            self.depth_optimization = False

        lb = np.append(lb, self.heating_section.optimization_bounds.lb)
        ub = np.append(ub, self.heating_section.optimization_bounds.ub)
        return Bounds(lb, ub)

    def set_optimization_param(self, optimization_param):

        if self.depth_optimization:

            self.dz_well = optimization_param[0]
            optimization_param = optimization_param[1:]

        self.heating_section.set_optimization_param(optimization_param)

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <---------------------   DISPLAY METHODS   ------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def __str__(self):

        class Formats:

            def __init__(self, n_element):

                table_width = 147
                col_width = int(table_width / n_element)
                actual_table_width = col_width * n_element

                self.__title = (ef.bold + "{:^" + str(actual_table_width) + "}" + ef.rs)
                self.__format_string = "{:^" + str(col_width) + "}"
                self.__format_first_column = ef.bold + "{:<" + str(col_width) + "}" + ef.rs
                self.__number_format_string = "{:.2f}"
                self.__header_format = ef.bold + self.__format_string + ef.rs
                self.__units_format = ef.italic + self.__format_string + ef.rs

            def format_number(self, number):

                try:

                    return self.__format_string.format(self.__number_format_string.format(number))

                except:

                    return self.__format_string.format(str(number))

            def format_header(self, header):
                return self.__header_format.format(header)

            def format_units(self, unit):
                return self.__units_format.format(unit)

            def format_title(self, title):
                return self.__title.format(title)

            def format_first_column(self, element):
                return self.__format_first_column.format(element)

        return "\n\n{}\n\n".format(

            self.__format_points_table(Formats)

        )

    def __format_points_table(self, formats):

        table_name = "POINTS"
        header_list = ["P", "T", "h", "s", "rho", "exergy"]

        frm = formats(len(header_list) + 1)

        rows_str = ""
        header_string = frm.format_first_column("Point")
        units_string = frm.format_first_column("-")
        counter = 0

        for point in self.points:

            row_str = frm.format_first_column(counter)
            counter += 1

            for element in header_list:

                if point == self.points[0]:

                    header_string += frm.format_header(element)
                    units_string += frm.format_units(point.get_unit(element))

                row_str += frm.format_number(point.get_variable(element))

            rows_str += "\n{}".format(row_str)

        return "\n{}\n\n{}\n{}\n{}".format(

            frm.format_title(table_name),
            header_string,
            units_string,
            rows_str

        )

    def __format_components_results(self, formats):

        table_name = "COMPONENTS EXERGY ANALYSIS"
        header_list = ["Name", "Exergy_fuel [kW]", "Exergy_product [kW]", "Exergy_destruction [kW]", "eta", "r", "y"]
        units_list = ["-", "kJ/kg", "kJ/kg", "kJ/kg", "-", "-", "-"]

        frm = formats(len(header_list))

        rows_str = ""
        header_string = ""
        units_string = ""

        for row in range(len(self.comp_results[header_list[0]])):

            row_str = ""
            for i in range(len(header_list)):

                element = header_list[i]

                if row == 0:

                    if "Exergy_" in element:

                        name = element.strip("Exergy_").strip(" [kW]")

                    else:

                        name = element

                    if not i == 0:

                        header_string += frm.format_header(name)
                        units_string += frm.format_units(units_list[i])

                    else:

                        header_string += frm.format_first_column(name)
                        units_string += frm.format_first_column(units_list[i])

                if i == 0:

                    row_str += frm.format_first_column(self.comp_results[element][row])

                else:

                    row_str += frm.format_number(self.comp_results[element][row])

            rows_str += "\n{}".format(row_str)

        return "\n{}\n\n{}\n{}\n{}".format(

            frm.format_title(table_name),
            header_string,
            units_string,
            rows_str

        )

    def get_iteration_profile(self, position_list):

        __tmp_point = self.points[0].duplicate()
        t_list = np.full((2, len(position_list)), np.nan)
        p_list = np.full((2, len(position_list)), np.nan)

        for i in range(len(position_list)):

            pos = position_list[i]

            p, rho = self.get_iteration_value(pos, 0)
            __tmp_point.set_variable("P", p)
            __tmp_point.set_variable("rho", rho)
            t_list[0, i] = __tmp_point.get_variable("T")
            p_list[0, i] = __tmp_point.get_variable("P")

            p, rho = self.get_iteration_value(self.dz_well - pos, 1)
            __tmp_point.set_variable("P", p)
            __tmp_point.set_variable("rho", rho)
            t_list[1, i] = __tmp_point.get_variable("T")
            p_list[1, i] = __tmp_point.get_variable("P")

        return np.array(t_list), np.array(p_list)

    def get_iteration_value(self, position, index):

        if len(self.integrators_profiler) > index:

            for step in self.integrators_profiler[index]:

                if step["range"][0] <= position <= step["range"][1]:

                    return step["dense_out"](position)
