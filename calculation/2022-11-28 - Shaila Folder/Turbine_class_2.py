# %% INIT
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
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

        self.evaluate_design()

    def evaluate_design(self):

        self.y_ids = list()
        self.dh_iso_des = list()
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
            self.points.append(self.input_point.duplicate())

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

            stage_output.set_variable("P", p_out)
            stage_output.set_variable("s", stage_input.get_variable("s"))

            h_iso_curr = stage_output.get_variable("h")
            delta_h_iso = h_in - h_iso_curr
            eta_curr = self.evaluate_eta(delta_h_iso, self.dh_iso_des[n])
            h_out = h_in - eta_curr * delta_h_iso

            stage_output.set_variable("P", p_out)
            stage_output.set_variable("h", h_out)

        if update_output:

            self.points[-1].copy_state_to(self.output_point)

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

    def evaluate_eta(self, delta_h_iso, delta_h_iso_des):

        a = np.log(delta_h_iso / delta_h_iso_des)
        return self.eta_des * (10 ** ((-0.00817 * (a ** 3)) - (0.03181 * (a ** 2)) + 0.0019 * a))

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

class TurbineODPlotter:

    def __init__(self, turbine: TurbineOD):

        self.turbine = turbine

    def plot_characteristic_curves(self, fig: Figure, n_points=20, off_design_conditions=None):

        ax_beta = fig.add_subplot(1, 2, 1)
        ax_eta = fig.add_subplot(1, 2, 2)

        fi_list, beta_list = self.__evaluate_beta_curve(n_points)
        dh_list, eta_list = self.__evaluate_eta_curve(n_points)

        ax_beta.plot(fi_list, beta_list)
        ax_eta.plot(dh_list, eta_list)

        if off_design_conditions is not None:

            dh_od, eta_od, fi_od, beta_od = self.__get_off_design_points(off_design_conditions)

            ax_beta.plot(fi_od, beta_od, "x")
            ax_eta.plot(dh_od, eta_od, "x")

    def __get_off_design_points(self, off_design_conditions):

        dh_od = list()
        eta_od = list()

        fi_od = list()
        beta_od = list()

        for off_design in off_design_conditions:

            self.turbine.update_input_output_points(off_design[0], off_design[1])
            self.turbine.update_off_design_flow_rate()

            fi, dh_iso, beta, eta = self.__get_od_plot_data()

            fi_od.append(fi)
            dh_od.append(dh_iso)
            beta_od.append(beta)
            eta_od.append(eta)

        return dh_od, eta_od, fi_od, beta_od

    def __evaluate_beta_curve(self, n_points=20):

        max_fi = self.turbine.get_max_fi()
        p_in = self.turbine.input_point.get_variable("P")
        fi_list = list()
        beta_list = list()

        for n in range(0, n_points):

            fi_perc = np.sqrt(float(n) / n_points)
            fi = max_fi * fi_perc

            self.turbine.input_point.m_dot = turbine.get_m_dot_from_fi(fi, self.turbine.input_point)
            self.turbine.update_off_design_output()

            if not math.isnan(self.turbine.points[-1].get_variable("P")):

                fi_list.append(fi)
                beta_list.append(p_in / self.turbine.output_point.get_variable("P"))

        return fi_list, beta_list

    def __evaluate_eta_curve(self, n_points=20):

        dh_iso_des = np.sum(self.turbine.dh_iso_des)
        dh_list = list()
        eta_list = list()

        for n in range(0, n_points):

            dh_perc = 1.2 * (float(n) + 1) / n_points + 0.4
            dh_real = dh_iso_des * dh_perc

            eta = self.turbine.evaluate_eta(dh_real, dh_iso_des)

            dh_list.append(dh_real)
            eta_list.append(eta)

        return dh_list, eta_list

    def __get_od_plot_data(self):

        p_in = self.turbine.input_point.get_variable("P")
        p_out = self.turbine.points[-1].get_variable("P")

        tmp_point = self.turbine.points[-1].duplicate()
        tmp_point.set_variable("P", p_out)
        tmp_point.set_variable("s", self.turbine.input_point.get_variable("s"))
        h_out_iso = tmp_point.get_variable("h")

        h_in = self.turbine.input_point.get_variable("h")
        h_out = self.turbine.points[-1].get_variable("h")
        rho_in =  self.turbine.input_point.get_variable("rho")

        beta = p_in / p_out
        dh_iso = h_in - h_out_iso
        eta = (h_in - h_out) / dh_iso
        fi = self.turbine.evaluate_fi(self.turbine.input_point.m_dot, p_in, rho_in)

        return fi, dh_iso, beta, eta


# %% DESIGN
from main_code.simplified_well.simplified_well import SimplifiedBHE

T_amb = 15 # [°C]
dT_appr = 7  # [°C]
delta_pump_p = 0.1

dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir


def evaluate_well():

    T_c_out=T_amb+dT_appr       # condenser Temperature
    CO2_input = PlantThermoPoint(["CarbonDioxide"], [1])
    CO2_input.set_variable("T", T_c_out)
    CO2_input.set_variable("Q", 0.)

    p_c = CO2_input.get_variable("P")
    CO2_input.set_variable("T", T_c_out)
    p_in_BH = p_c + delta_pump_p
    CO2_input.set_variable("P", p_in_BH)

    bhe_inside = SimplifiedBHE(

        input_thermo_point=CO2_input,
        dz_well=dz_well, T_rocks=T_rock, use_rk=True

    )

    bhe_inside.update()

    return bhe_inside.points[-1].duplicate(), CO2_input.duplicate()

turbine_in, turbine_out = evaluate_well()
turbine = TurbineOD(turbine_in, turbine_out, n_stages=3)

# %% OFF-DESIGN
T_amb = 15 # [°C]
dT_appr = 7  # [°C]

beta_OD_list=list()
Fi_OD_list=list()
Power_OD_list=list()
m_dot_list=list()
T_a=list()
for T in range(0,24):
    T_amb=T
    T_a.append(T)
    turbine_in, turbine_out = evaluate_well()
    turbine.update_input_output_points(turbine_in, turbine_out)
    turbine.update_off_design_flow_rate()

    m_dot_od = turbine.input_point.m_dot
    m_dot_list.append(m_dot_od)

    h_in_od=turbine_in.get_variable("h")
    h_out_od=turbine_out.get_variable("h")
    rho_in=turbine_in.get_variable("rho")

    P_in=turbine_in.get_variable("P")
    P_out=turbine_out.get_variable("P")
    b=P_in/P_out
    Pow=(h_in_od-h_out_od)*m_dot_od
    Power_OD_list.append(Pow)
    beta_OD_list.append(b)
    fi=m_dot_od/(np.sqrt(P_in*rho_in))
    Fi_OD_list.append(fi)

# %% PLOT
plt.plot(T_a,Fi_OD_list)
plt.xlabel('T_amb/°c')
plt.ylabel('Mass Flow')
plt.show()
# %% PLOT
plt.plot(T_a,m_dot_list,'r')
plt.xlabel('T_amb/°c')
plt.ylabel('m_dot')
plt.show()
# %% PLOT
plt.plot(T_a,Power_OD_list,'g')
plt.xlabel('T_amb/°c')
plt.ylabel('Power/W')
plt.show()
# %% PLOT
plt.plot(T_a,beta_OD_list,'k')
plt.xlabel('T_amb/°c')
plt.ylabel('beta')
plt.show()

# %% INPUT POINTS

conditions = list()
for T in range(-5, 20):

    T_amb = T
    turbine_in, turbine_out = evaluate_well()
    conditions.append([turbine_in, turbine_out])


# %% PLOT

fig = plt.figure(dpi=150)
fig.set_size_inches(20, 8)
turbine_plotter = TurbineODPlotter(turbine)
turbine_plotter.plot_characteristic_curves(fig, n_points=30, off_design_conditions=conditions)
plt.tight_layout()
plt.show()
# %% PLOT
fig.set_size_inches(20, 8)
turbine_plotter = TurbineODPlotter(turbine)
