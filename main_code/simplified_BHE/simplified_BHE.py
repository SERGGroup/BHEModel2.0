from main_code.simplified_BHE.heating_sections import DefaultHeatingSection, AbstractHeatingSection
from main_code.support import PlantThermoPoint, retrieve_PPI
from main_code import constants

from EEETools.Tools.API.ExcelAPI.modules_importer import calculate_excel
from EEETools.Tools.API.Tools.main_tools import get_result_data_frames
from scipy.optimize import Bounds
from scipy.constants import g
import numpy as np
from sty import ef
import warnings
import os


class SimplifiedBHE:

    def __init__(

            self, input_thermo_point: PlantThermoPoint, dz_well, T_rocks,
            heating_section=None, k_rocks=0.2, geo_flux=0.1, q_up_down=0.0,
            PPI=None

    ):

        # GEOTHERMAL FIELD CONDITION

        self.dz_well = dz_well      # unit: m .............. default: NONE (REQUIRED)
        self.T_rocks = T_rocks      # unit: °C ............. default: NONE (REQUIRED)
        self.k_rocks = k_rocks      # unit: W / (m K) ...... default: 0.2 W / (m K)
        self.geo_flux = geo_flux    # unit: W / m^2 ........ default: 0.1 W / m^2

        self.geo_gradient = None
        self.depth_optimization = False

        self.q_dot_up = - q_up_down
        self.q_dot_down = q_up_down

        self.c_well = 0.
        self.C0 = [0., 0.]

        self.__init_points(input_thermo_point)
        self.__init_heating_section(heating_section)
        self.__reset_control_elements(first_initialization=True)
        self.__ambient_condition = None

        self.__current_PPI = PPI

    def __init_points(self, input_thermo_point):

        self.points = [input_thermo_point]

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

        self.__evaluate_points()
        self.__evaluate_parameters()
        self.__evaluate_performance_parameters()
        self.__perform_exergy_analysis()

        self.__reset_control_elements()

    def __set_reference_state(self, ambient_condition):

        self.__ambient_condition = ambient_condition

        if type(ambient_condition) == PlantThermoPoint:

            ref_state = ambient_condition

        else:

            ref_state = self.points[0].duplicate()

        self.points[0].reference_state = ref_state

    def __evaluate_points(self):

        self.C0[0] = self.__update_DP_vertical(self.points[0], is_upward=False)
        self.heating_section.update()
        self.C0[1] = self.__update_DP_vertical(self.points[2], is_upward=True)

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

    def __update_DP_vertical(self, input_point: PlantThermoPoint, is_upward=True):

        if input_point == self.points[-1]:

            new_point = input_point.duplicate()
            self.points.append(new_point)

        else:

            new_point = self.points[self.points.index(input_point) + 1]

        if is_upward:

            mult = -1
            dq = self.q_dot_up

        else:

            mult = 1
            dq = self.q_dot_down

        rho_in = input_point.get_variable("rho")
        dP = mult * (rho_in * self.dz_well * g) / 10 ** 6
        Dh = (mult * g / 10 ** 3 + dq) * self.dz_well
        dh_dp_stream = (g + dq * 10 ** 3) / g * 10 ** 3

        new_point.set_variable("h", input_point.get_variable("h") + Dh)
        new_point.set_variable("P", input_point.get_variable("P") + dP)

        CO_in = self.__calculate_C0(input_point, dh_dp_stream)
        CO_out = self.__calculate_C0(new_point, dh_dp_stream)
        CO = (CO_out + CO_in) / 2
        counter = 0

        while True:

            dP = mult * abs((1 - np.exp(CO * g * self.dz_well / 10 ** 6)) * rho_in / CO)
            new_point.set_variable("P", input_point.get_variable("P") + dP)

            old_CO = CO

            CO_out = self.__calculate_C0(new_point, dh_dp_stream)
            CO = (CO_out + CO_in) / 2

            if abs((old_CO - CO) / CO) < 10 ** -10 or counter > 20:

                break

            else:

                counter += 1

        return CO

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