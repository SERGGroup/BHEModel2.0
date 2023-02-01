# %% INIT
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import numpy as np
import math

class TurbineOD:

    def __init__(

            self,
            input_point:PlantThermoPoint, output_point:PlantThermoPoint,
            n_stages=2, eta_des=0.7

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
            fi = self.input_point.m_dot / np.sqrt(p_in * rho_in)

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

            fi_stage = m_dot / np.sqrt(p_in * rho_in)
            beta_stage = 1 / np.sqrt(1 - ((fi_stage ** 2) * self.y_ids[n]))
            p_out = p_in / beta_stage

            stage_output.set_variable("P", p_out)
            stage_output.set_variable("s", stage_input.get_variable("s"))

            h_iso_curr = stage_output.get_variable("h")
            delta_h_iso = h_in - h_iso_curr
            a = np.log(delta_h_iso / self.dh_iso_des[n])
            eta_curr = self.eta_des * (10 ** ((-0.00817 * (a ** 3)) - (0.03181 * (a ** 2)) + 0.0019 * a))
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

        while counter < 20:

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

    def update_input_output_points(self, turbine_in, turbine_out):

        turbine_in.copy_state_to(self.input_point)
        turbine_out.copy_state_to(self.output_point)


# %% DESIGN
from main_code.simplified_well.simplified_well import SimplifiedBHE

T_amb = 15 # [°C]
dT_appr = 7  # [°C]
delta_pump_p = 0.1

dz_well = 1500  # [m] Depth of the reservoir
T_rock = 125    # [°C] Temperature of the reservoir

T_c_out=T_amb+dT_appr       # condenser Temperature

def evaluate_well():

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

    return bhe_inside.points[-1], CO2_input

turbine_in, turbine_out = evaluate_well()
turbine = TurbineOD(turbine_in, turbine_out, n_stages=3)

# %% OFF-DESIGN
T_amb = 15 # [°C]
dT_appr = 7  # [°C]

turbine_in, turbine_out = evaluate_well()
turbine.update_input_output_points(turbine_in, turbine_out)
turbine.update_off_design_flow_rate()
print(turbine.input_point.m_dot)