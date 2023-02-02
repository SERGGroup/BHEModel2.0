from main_code.off_design_model.turbine_off_design import TurbineOD
from matplotlib.figure import Figure
import numpy as np
import math

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

        ax_beta.set_xlabel(r'$\phi$ [-]', fontsize='large', loc='center')
        ax_beta.set_ylabel(r'$\beta$ [-]', fontsize='large', loc='center')
        ax_beta.set_title("Stodola Curve\n$^{\it{(overall)}}$", fontsize='xx-large', loc='center')

        ax_eta.set_xlabel(r'$dh_{iso}$ [kJ/kg]', fontsize='large', loc='center')
        ax_eta.set_ylabel(r'$\eta_{iso}$ [-]', fontsize='large', loc='center')
        ax_eta.set_title("Efficiency Curve\n$^{\it{(1^{st}\ stage)}}$", fontsize='xx-large', loc='center')

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

            self.turbine.input_point.m_dot = self.turbine.get_m_dot_from_fi(fi, self.turbine.input_point)
            self.turbine.update_off_design_output()

            if not math.isnan(self.turbine.points[-1].get_variable("P")):

                fi_list.append(fi)
                beta_list.append(p_in / self.turbine.output_point.get_variable("P"))

        return fi_list, beta_list

    def __evaluate_eta_curve(self, n_points=20):

        dh_iso_des = self.turbine.dh_iso_des[0]
        dh_list = list()
        eta_list = list()

        for n in range(0, n_points):

            dh_perc = 1.2 * (float(n) + 1) / n_points + 0.4
            dh_real = dh_iso_des * dh_perc

            eta = self.turbine.evaluate_eta_stage(dh_real, dh_iso_des)

            dh_list.append(dh_real)
            eta_list.append(eta)

        return dh_list, eta_list

    def __get_od_plot_data(self):

        p_in = self.turbine.input_point.get_variable("P")
        p_out = self.turbine.points[-1].get_variable("P")

        tmp_point = self.turbine.points[-1].duplicate()
        tmp_point.set_variable("P", self.turbine.points[1].get_variable("P"))
        tmp_point.set_variable("s", self.turbine.input_point.get_variable("s"))
        h_1st_out_iso = tmp_point.get_variable("h")

        h_in = self.turbine.input_point.get_variable("h")
        h_1st_out = self.turbine.points[1].get_variable("h")
        rho_in =  self.turbine.input_point.get_variable("rho")

        beta = p_in / p_out
        dh_1st_iso = h_in - h_1st_out_iso
        eta_1st = (h_in - h_1st_out) / dh_1st_iso
        fi = self.turbine.evaluate_fi(self.turbine.input_point.m_dot, p_in, rho_in)

        return fi, dh_1st_iso, beta, eta_1st
