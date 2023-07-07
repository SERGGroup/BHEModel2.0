from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from sty import ef
import numpy as np
import warnings
import math

class TurbineOD:

    def __init__(

            self,
            input_point:PlantThermoPoint, output_point:PlantThermoPoint,
            n_stages=3, eta_des=0.7

    ):

        self.input_point = input_point
        self.output_point = output_point
        self.n_stages = n_stages
        self.eta_des = eta_des

        self.y_ids = list()
        self.dh_iso_des = list()

        self.points = list()
        self.design_points = list()

        self.power = 0.
        self.eta_iso = 0.

        self.evaluate_design()

    def evaluate_design(self):

        self.y_ids = list()
        self.dh_iso_des = list()
        self.lf_des = list()

        self.points = [self.input_point]
        self.design_points = [self.input_point.duplicate()]

        p_levels = self.__evaluate_p_levels()

        for n in range(self.n_stages):

            stage_input = self.design_points[-1]

            stage_output = stage_input.duplicate()
            stage_output.set_variable("P", p_levels[n])
            stage_output.set_variable("s", stage_input.get_variable("s"))

            dh_is = stage_input.get_variable("h") - stage_output.get_variable("h")
            h_out = stage_input.get_variable("h") - dh_is * self.eta_des

            stage_output.set_variable("P", p_levels[n])
            stage_output.set_variable("h", h_out)

            p_in = stage_input.get_variable("P")
            rho_in = stage_input.get_variable("rho")
            beta = stage_input.get_variable("P") / stage_output.get_variable("P")
            fi = self.evaluate_fi(self.input_point.m_dot, p_in, rho_in)

            self.dh_iso_des.append(dh_is)
            self.y_ids.append((1 - (1 / (beta ** 2))) / fi ** 2)

            self.design_points.append(stage_output)
            self.points.append(stage_output.duplicate())

        self.__evaluate_parameters()

    def __evaluate_p_levels(self):

        tmp_point = self.input_point.duplicate()
        tmp_point.set_variable("P", self.output_point.get_variable("P"))
        tmp_point.set_variable("s", self.input_point.get_variable("s"))

        s_t_in = self.input_point.get_variable("s")
        h_t_in = self.input_point.get_variable("h")
        h_out_iso = tmp_point.get_variable("h")

        p_levels = list()
        delta_h_t_iso = h_t_in - h_out_iso
        delta_h_t_i = delta_h_t_iso / self.n_stages

        for i in range(self.n_stages):

            tmp_point.set_variable("h", h_t_in - delta_h_t_i * (i + 1))
            tmp_point.set_variable("s", s_t_in)
            p_levels.append(tmp_point.get_variable("P"))

        return p_levels

    def __evaluate_parameters(self):

        tmp_point = self.points[-1].duplicate()
        tmp_point.set_variable("P", self.points[-1].get_variable("P"))
        tmp_point.set_variable("s", self.points[0].get_variable("s"))

        dh = self.points[0].get_variable("h") - self.points[-1].get_variable("h")
        dh_iso = self.points[0].get_variable("h") - tmp_point.get_variable("h")

        self.power = dh * self.input_point.m_dot
        self.eta_iso = dh / dh_iso

    def update_off_design_output(self, update_output=True):

        for n in range(len(self.y_ids)):

            stage_input = self.points[n]
            stage_output = self.points[n + 1]
            m_dot = self.input_point.m_dot
            stage_output.m_dot = m_dot

            p_in = stage_input.get_variable("P")
            h_in = stage_input.get_variable("h")
            rho_in = stage_input.get_variable("rho")

            fi_stage = self.evaluate_fi(m_dot, p_in, rho_in)
            beta_stage = 1 / np.sqrt(1 - ((fi_stage ** 2) * self.y_ids[n]))
            p_out = p_in / beta_stage

            eta_curr, delta_h_iso = self.evaluate_eta_stage(n, stage_input, p_out)
            h_out = h_in - eta_curr * delta_h_iso

            stage_output.set_variable("P", p_out)
            stage_output.set_variable("h", h_out)

        if update_output:

            self.points[-1].copy_state_to(self.output_point)
            self.__evaluate_parameters()

    def update_off_design_flow_rate(self):

        counter = 0
        m_r_low = 0
        m_r_high = 2.5
        p_out_target = self.output_point.get_variable("p")

        warnings.filterwarnings('ignore')  # ignore warning

        while counter < 40:

            m_r = (m_r_high + m_r_low) / 2
            m_dot_guess = m_r * self.design_points[0].m_dot

            self.input_point.m_dot = m_dot_guess
            self.update_off_design_output(update_output=False)

            if math.isnan(self.points[-1].get_variable("P")):

                m_r_high = m_r

            else:

                if p_out_target - self.points[-1].get_variable("P") < 0:

                    m_r_low = m_r

                else:

                    m_r_high = m_r

            counter += 1

        warnings.filterwarnings('default')  # restore warning

        m_r = m_r_low
        m_dot_guess = m_r * self.design_points[0].m_dot

        self.input_point.m_dot = m_dot_guess
        self.update_off_design_output(update_output=False)
        self.__evaluate_parameters()

    def update_input_output_points(self, turbine_in, turbine_out):

        turbine_in.copy_state_to(self.input_point)
        turbine_out.copy_state_to(self.output_point)

    def get_max_fi(self):

        counter = 0
        m_r_low = 0
        m_r_high = 2.5
        m_r = (m_r_high + m_r_low) / 2

        while counter < 20:

            m_r = (m_r_high + m_r_low) / 2
            m_dot_guess = m_r * self.design_points[0].m_dot
            self.input_point.m_dot = m_dot_guess
            self.update_off_design_output(update_output=False)

            if math.isnan(self.points[-1].get_variable("P")):

                m_r_high = m_r

            else:

                m_r_low = m_r

            counter += 1

        return self.evaluate_fi(m_r, self.input_point.get_variable("P"), self.input_point.get_variable("rho"))

    def evaluate_eta_stage(self, n, stage_input, p_out):

        tmp_point= stage_input.duplicate()
        tmp_point.set_variable("P", p_out)
        tmp_point.set_variable("s", stage_input.get_variable("s"))

        h_iso_curr = tmp_point.get_variable("h")
        dh_iso = stage_input.get_variable("h") - h_iso_curr

        a = np.log(dh_iso / self.dh_iso_des[n])
        eta_stage = self.eta_des * (10 ** ((-0.00817 * (a ** 3)) - (0.03181 * (a ** 2)) + 0.0019 * a))

        return eta_stage, dh_iso

    @staticmethod
    def evaluate_fi(m_dot, p_in, rho_in):

        return m_dot / np.sqrt(p_in * 1E6 * rho_in)

    @staticmethod
    def get_m_dot_from_fi(fi, point:PlantThermoPoint):

        p_in = point.get_variable("P") * 1E6
        rho_in = point.get_variable("rho")
        return fi * np.sqrt(p_in * rho_in)

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