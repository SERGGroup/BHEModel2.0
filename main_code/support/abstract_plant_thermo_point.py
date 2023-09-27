from REFPROPConnector import ThermodynamicPoint
from sty import ef


class PlantThermoPoint(ThermodynamicPoint):

    def __init__(

            self, fluids: list, composition: list, m_dot=1., rp_handler=None,
            other_variables="all", calculate_on_need="all", unit_system="SI WITH C"

    ):

        super().__init__(fluids, composition, rp_handler, other_variables, calculate_on_need, unit_system)
        self.m_dot = m_dot

    def duplicate(self):

        tp = PlantThermoPoint(

            self.RPHandler.fluids,
            self.RPHandler.composition,
            self.m_dot,
            rp_handler=self.RPHandler,
            unit_system=self.RPHandler.unit_system,
            other_variables=self.inputs["other_variables"],
            calculate_on_need=self.inputs["calculate_on_need"]

        )

        self.copy_state_to(tp)
        return tp

    def copy_state_to(self, target_point):

        super().copy_state_to(target_point)
        target_point.m_dot = self.m_dot

    def set_to_compression_result(self, P_out, eta, input_state: ThermodynamicPoint):

        is_expansion = input_state.get_variable("P") > P_out

        self.set_variable("P", P_out)
        self.set_variable("s", input_state.get_variable("s"))

        h_iso = self.get_variable("h")
        h_in = input_state.get_variable("h")

        if is_expansion:

            h_out = h_in - (h_in - h_iso) * eta

        else:

            h_out = (h_iso - h_in) / eta + h_in

        self.set_variable("h", h_out)
        self.set_variable("P", P_out)

    def set_to_expansion_result(self, P_out, eta, input_state: ThermodynamicPoint):

        self.set_to_compression_result(P_out, eta, input_state)

    def set_to_specific_super_heating(self, P_value, super_heating_value):

        self.set_variable("P", P_value)
        self.set_variable("x", 0)

        T_sat = self.get_variable("T")
        T_new = T_sat + super_heating_value

        self.set_variable("P", P_value)
        self.set_variable("T", T_new)

    def dvar(self, input_state: ThermodynamicPoint, variable):

        """
        Evaluate Variable "variable" variation with respect to "input_state".
        For example: point.dvar(other_point, "h") will evaluate the enthalpy variation between point and other_point
        """
        return self.get_variable(variable) - input_state.get_variable(variable)

    def dex(self, input_state: ThermodynamicPoint, T_ref):
        """
        Evaluate exergy variation with respect to "input_state". using "T_ref" as the reference temperature
        """
        dh = self.dvar(input_state, "h")
        ds = self.dvar(input_state, "s")

        if self.get_unit("T") == "K":

            return dh - T_ref * ds

        else:

            return dh - (T_ref + 273.15) * ds

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

        table_name = self.str_fluid
        header_list = ["P", "T", "h", "s", "rho", "Q"]

        frm = formats(len(header_list) + 1)

        row_str = ""
        header_string = ""
        units_string = ""

        for element in header_list:

            header_string += frm.format_header(element)
            units_string += frm.format_units(self.get_unit(element))

            if not element == "Q":

                row_str += frm.format_number(self.get_variable(element))

            else:

                Q_value = self.get_variable(element)

                if Q_value is not None and 0 <= Q_value <= 1:

                    row_str += frm.format_number(self.get_variable(element))

                else:

                    row_str += frm.format_units(" - ")

        return "\n{}\n\n{}\n{}\n{}".format(

            frm.format_title(table_name),
            header_string,
            units_string,
            row_str

        )

    @property
    def str_fluid(self):

        fluid_list = self.RPHandler.fluids
        comp_list = self.RPHandler.composition

        if len(fluid_list) == 1:

            return fluid_list[0]

        else:

            form_str = "{:.0f}% {}"
            fluid_name = form_str.format(comp_list[0] * 100, fluid_list[0])

            for i in range(1, len(fluid_list)):

                fluid_name += " - {}".format(form_str.format(comp_list[i] * 100, fluid_list[i]))

            return fluid_name