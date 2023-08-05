from main_code.well_model.simplified_well.heating_sections import DefaultHeatingSection, AbstractHeatingSection
from main_code.support.other.integration_profiler import IntegrationProfiler
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

            self, input_thermo_point: PlantThermoPoint, dz_well,
            t_rocks=None, t_surf=None, geo_flux=None,
            k_rocks=0.2, c_rocks=1, rho_rocks=2500,
            heating_section=None, PPI=None,
            use_rk=True

    ):

        # GEOTHERMAL FIELD CONDITION
        self.__init_geological_data(dz_well, t_rocks, t_surf, geo_flux, k_rocks, c_rocks, rho_rocks)

        self.depth_optimization = False
        self.use_rk = use_rk

        self.c_well = 0.
        self.C0 = [0., 0.]

        self.integrators_profiler = list()
        self.__init_points(input_thermo_point)
        self.__init_heating_section(heating_section)
        self.__reset_control_elements(first_initialization=True)

        self.__ambient_condition = None
        self.is_upward = False

        self.__current_PPI = PPI

    def __init_geological_data(self, dz_well, t_rocks, t_surf, geo_flux, k_rocks, c_rocks, rho_rocks):

        self.dz_well = dz_well  # unit: m .............. default: NONE (REQUIRED)
        self.k_rocks = k_rocks  # unit: W / (m K) ...... default: 0.2 W / (m K)
        self.c_rocks = c_rocks  # unit: kJ / (kg K) .... default: 1 kJ / (kg K)
        self.rho_rocks = rho_rocks  # unit: kg / m^3 ....... default: 2500 kg / m^3

        #   alpha_rocks -> rock thermal diffusivity in [m^2/s]
        #   (1e3 conversion factor c_rocks [kJ / (kg K)] -> [J / (kg K)])
        self.alpha_rocks = self.k_rocks / (self.rho_rocks * self.c_rocks * 1e3)
        self.__init_temperatures(t_rocks, t_surf, geo_flux)

    def __init_temperatures(self, t_rocks, t_surf, geo_flux):

        if geo_flux is None or (t_surf is not None and t_rocks is not None):

            if t_surf is None:

                t_surf = 10 # [°C] (default value)

            if t_rocks is None:

                self.__raise_t_rocks_error()

            self.t_surf = t_surf
            self.t_rocks = t_rocks
            self.geo_gradient = (self.t_rocks - self.t_surf) / self.dz_well
            self.geo_flux = self.geo_gradient * self.k_rocks

        elif t_surf is None:

            if geo_flux is None:

                geo_flux = 0.1

            if t_rocks is None:

                self.__raise_t_rocks_error()

            self.t_rocks = t_rocks
            self.geo_flux = geo_flux    # unit: W / m^2 ........ default: 0.1 W / m^2
            self.geo_gradient = self.geo_flux / self.k_rocks
            self.t_surf = self.t_rocks - self.geo_gradient * self.dz_well

        elif t_rocks is None:

            if geo_flux is None or t_surf is None:

                self.__raise_t_rocks_error()

            self.t_surf = t_surf
            self.geo_flux = geo_flux  # unit: W / m^2 ........ default: 0.1 W / m^2
            self.geo_gradient = self.geo_flux / self.k_rocks
            self.t_rocks = self.t_surf + self.geo_gradient * self.dz_well

    def __init_points(self, input_thermo_point):

        self.points = [input_thermo_point]
        self.new_point = input_thermo_point.duplicate()

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

    def __raise_t_rocks_error(self):

        raise ValueError(

            """

                \n't_rock' is needed to start the calculations: 
                \nAs an alternative, you can provide both 't_surf' and 'geo_flux'

            """

        )

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

    @property
    def input_point(self):

        return self.points[0]

    @property
    def output_point(self):

        return self.points[-1]

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

        if self.q_bottom > 0:
            self.eta_I = self.dh / self.q_bottom

        else:
            self.eta_I = 0

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

        if self.q_bottom > 0:

            excel_path = os.path.join(constants.RES_FOLDER, "3ETool_Excel-Files", "BHE_simple_analysis.xlsx")
            T_rocks_K = self.t_rocks + 273.15

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

        else:

            self.comp_results = 0.
            self.eta_II = 0.

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------   VERTICAL WELL INTEGRATION   ------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update_DP_vertical(self, input_point: PlantThermoPoint, is_upward=True):

        if input_point == self.points[-1]:
            self.new_point = input_point.duplicate()
            self.points.append(self.new_point)

        else:
            self.new_point = self.points[self.points.index(input_point) + 1]

        self.is_upward = is_upward
        if is_upward:
            depths = (self.dz_well, 0)

        else:
            depths = (0, self.dz_well)

        p0 = input_point.get_variable("P")
        rho0 = input_point.get_variable("rho")

        if not self.use_rk:

            integrator = IntegrationProfiler(

                SimpleIntegrator(self.rk_funct, depths[0], [p0, rho0], depths[1], n_steps=100)

            )

        else:

            integrator = IntegrationProfiler(

                RK45(self.rk_funct, depths[0], [p0, rho0], depths[1])

            )

        self.integrators_profiler.append(integrator)

        while integrator.status == 'running':

            integrator.step()

        output = integrator.y
        self.new_point.set_variable("P", output[0])
        self.new_point.set_variable("rho", output[1])

        co_in = self.calculate_C0(input_point, depths[0])
        co_out = self.calculate_C0(self.new_point, depths[1])
        c0 = (co_in + co_out) / 2

        self.additional_final_calculations(input_point, depths)
        self.new_point = self.new_point.duplicate()

        return c0

    def rk_funct(self, z, y):

        p_curr = y[0]
        rho_curr = y[1]

        self.new_point.set_variable("P", p_curr)
        self.new_point.set_variable("rho", rho_curr)

        dp_loss = self.__evaluate_pressure_losses(self.new_point)
        dp = rho_curr * g / 1e6 + dp_loss

        if self.is_upward:
            c0_curr = self.calculate_C0(self.new_point, depth=self.dz_well - z, dp_dl=dp)

        else:
            c0_curr = self.calculate_C0(self.new_point, depth=z, dp_dl=dp)

        d_rho = c0_curr * dp

        return [dp, d_rho]

    def __evaluate_pressure_losses(self, input_point):

        dp = self.evaluate_pressure_losses(input_point)

        if self.is_upward:

            return dp

        else:

            return -dp

    def calculate_C0(self, curr_point: PlantThermoPoint, depth=0., dp_dl=None):

        drho_dP = curr_point.get_derivative("rho", "P", "T")
        drho_dT = curr_point.get_derivative("rho", "T", "P")
        dh_dT = curr_point.get_derivative("h", "T", "P")
        dh_dP = curr_point.get_derivative("h", "P", "T")

        if dp_dl is None:
            dp_dl = curr_point.get_variable("rho") * g / 1e6

        return drho_dP + drho_dT / dh_dT * (self.dh_dl_stream(curr_point, depth) / dp_dl - dh_dP)

    # The following three methods can be overwritten in subclasses to implement pressure
    # and heat transfer calculations

    def evaluate_pressure_losses(self, curr_point: PlantThermoPoint):

        return 0.

    def dh_dl_stream(self, curr_point: PlantThermoPoint, depth=0.):

        return g / 1e3

    def additional_final_calculations(self, input_point, depths):

        pass

    @property
    @abstractmethod
    def neglect_internal_heat_transfer(self):

        return False

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

    def get_iteration_points(self, position_list, profile=None):

        if profile is None:
            profile = self.integrators_profiler

        __tmp_point = self.points[0].duplicate()
        point_in_list = list()
        point_out_list = list()

        for i in range(len(position_list)):

            pos = position_list[i]

            p, rho = profile[0].get_iteration_value(pos)
            __tmp_point.set_variable("P", p)
            __tmp_point.set_variable("rho", rho)
            point_in_list.append(__tmp_point.duplicate())

            p, rho = profile[1].get_iteration_value(pos)
            __tmp_point.set_variable("P", p)
            __tmp_point.set_variable("rho", rho)
            point_out_list.append(__tmp_point.duplicate())

        return [point_in_list, point_out_list]

    def get_iteration_profile(self, position_list, profile=None):

        if profile is None:
            profile = self.integrators_profiler

        __tmp_point = self.points[0].duplicate()
        t_list = np.full((2, len(position_list)), np.nan)
        p_list = np.full((2, len(position_list)), np.nan)

        for i in range(len(position_list)):

            pos = position_list[i]

            p, rho = profile[0].get_iteration_value(pos)
            __tmp_point.set_variable("P", p)
            __tmp_point.set_variable("rho", rho)
            t_list[0, i] = __tmp_point.get_variable("T")
            p_list[0, i] = __tmp_point.get_variable("P")

            p, rho = profile[1].get_iteration_value(pos)
            __tmp_point.set_variable("P", p)
            __tmp_point.set_variable("rho", rho)
            t_list[1, i] = __tmp_point.get_variable("T")
            p_list[1, i] = __tmp_point.get_variable("P")

        return np.array(t_list), np.array(p_list)

    @property
    def calculation_setup_data(self):

        data_frame = {

            "Reservoir Data": {

                "depth": {"value": self.dz_well, "unit": "m"},
                "T_rock": {"value": self.t_rocks, "unit": "°C"},
                "T_surf": {"value": self.t_surf, "unit": "°C"},
                "geo_flux": {"value": self.geo_flux, "unit": "W / m^2"},
                "geo_gradient": {"value": self.geo_gradient, "unit": "°C/m"},

                "k_rocks": {"value": self.rho_rocks, "unit": "W / (m K)"},
                "c_rocks": {"value": self.c_rocks, "unit": "kJ / (kg K)"},
                "rho_rocks": {"value": self.rho_rocks, "unit": "kg / m^3"},
                "alpha_rocks": {"value": self.alpha_rocks, "unit": "-"}

            },

            "Input Conditions": {

                "fluid": {"value": self.points[0].str_fluid, "unit": None},
                "P_in": {"value": self.points[0].get_variable("P"), "unit": self.points[0].get_unit("P")},
                "T_in": {"value": self.points[0].get_variable("T"), "unit": self.points[0].get_unit("T")},
                "rho_in": {"value": self.points[0].get_variable("rho"), "unit": self.points[0].get_unit("rho")},
                "m_dot": {"value": self.points[0].m_dot, "unit": "kg/s"},

            },

            "Calculation Options": {

                "use rk": {"value": self.use_rk, "unit": None},
                "well class": {"value": "SimplifiedWell", "unit": None},
                "well model": {"value": "None", "unit": None},
                "pressure losses": {"value": "ignored", "unit": None},
                "heat transfer": {"value": "ignored", "unit": None}

            },

            "well geometry": {

                "depth": {"value": self.dz_well, "unit": "m"},

            },

            "Heating Section Data": {


            }

        }

        data_frame = self.additional_setup_data(data_frame)
        data_frame = self.heating_section.additional_setup_data(data_frame)
        return data_frame

    @abstractmethod
    def additional_setup_data(self, data_frame: dict):

        return data_frame


class SimplifiedBHE(SimplifiedWell):

    def evaluate_points(self):

        self.C0[0] = self.update_DP_vertical(self.points[0], is_upward=False)
        self.heating_section.update()
        self.C0[1] = self.update_DP_vertical(self.points[2], is_upward=True)

    def additional_setup_data(self, data_frame: dict):

        data_frame["Calculation Options"].update({

            "well class": {"value": "SimplifiedBHE", "unit": None},
            "well model": {"value": "BHE", "unit": None},

        })
        return data_frame

    @property
    def neglect_internal_heat_transfer(self):
        return False


class SimplifiedCPG(SimplifiedWell):

    def __init__(

            self, input_thermo_point: PlantThermoPoint, dz_well,
            t_rocks=None, t_surf=None, geo_flux=None,
            k_rocks=0.2, c_rocks=1, rho_rocks=2500,
            heating_section=None, PPI=None,
            use_rk=True, p_res=None

    ):

        super().__init__(

            input_thermo_point, dz_well, t_rocks=t_rocks, t_surf=t_surf,
            heating_section=heating_section, k_rocks=k_rocks, c_rocks=c_rocks,
            rho_rocks=rho_rocks, geo_flux=geo_flux, PPI=PPI, use_rk=use_rk

        )

        if p_res is not None:

            self.p_res = p_res

        else:

            self.p_res = 1000 * 9.81 * dz_well / 1e6  # water hydrostatic pressure (in MPa)

    def evaluate_points(self):

        # iterate point 0 in order to get fixed point 1 (production well bottom hole) pressure.
        # pressure = reservoir pressure - pressure losses in the reservoir

        counter = 0
        surface_temperature = self.points[0].get_variable("T")

        while True:

            self.C0[0] = self.update_DP_vertical(self.points[0], is_upward=False)
            self.heating_section.update()
            dp = self.points[2].get_variable("P") - self.p_res

            if abs(dp / self.p_res) < 1e-3 or counter > 10:

                break

            else:

                counter += 1
                self.points[0].set_variable("P", self.points[0].get_variable("P") - dp)
                self.points[0].set_variable("T", surface_temperature)

        self.C0[1] = self.update_DP_vertical(self.points[2], is_upward=True)

    def additional_setup_data(self, data_frame: dict):

        data_frame["Calculation Options"].update({

            "well class": {"value": "SimplifiedCPG", "unit": None},
            "well model": {"value": "CPG", "unit": None},

        })

        data_frame["Reservoir Data"].update({

            "P_reservoir": {"value": self.p_res, "unit": "MPa"}

        })

        return data_frame

    @property
    def neglect_internal_heat_transfer(self):
        return False