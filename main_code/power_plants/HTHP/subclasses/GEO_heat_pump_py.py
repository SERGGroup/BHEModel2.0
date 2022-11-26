from main_code.power_plants.HTHP.abstract_hthp import AbstractPythonHTHP
from main_code.power_plants.support import TURTON_PEC_calculator
from main_code.support import PlantThermoPoint
from scipy.optimize import Bounds
import numpy as np


class DirectWaterHeatPumpThermo(AbstractPythonHTHP):

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def optimization_bounds(self):

        return Bounds([0.4], [0.9])

    @property
    def initial_guess(self):

        return [0.6]

    def set_optimization_param(self, optimization_param):

        self.T_sc_perc = optimization_param[0]

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   INITIALIZATION METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def initialize_other_parameter(self):

        self.__T_sc_perc = 0.4  # Value to be optimized that control the water pressure (Limited Between 0-1)
        self.__n_stages_comp = 3  # Inter-cooled Compression stages
        self.IC_dt_SH = 10  # Residual super-heating after inter-cooling

        self.__init_fixed_parameters()
        self.__evaluate_T_sc()
        self.__init_points()

    def __evaluate_T_sc(self):

        T_sc_max = self.BHE_output.get_variable("T")
        T_sc_min = 10  # °C

        self.T_sc = (T_sc_max - T_sc_min) * self.T_sc_perc + T_sc_min

    def __init_points(self):

        self.points = list()

        # Point 0 - BHE Output

        self.points.append(self.BHE_output)

        # Points 1-7 - steam distribution section:
        #
        #   1 - Steam Cylinder Inlet
        #   2 - Steam Cylinder Vapour Outlet
        #   3 - Inter-cooled Compression Section Outlet
        #   4 - Steam Send
        #   5 - Steam Return
        #   6 - RH outlet
        #   7 - Steam Cylinder Steam Return

        for i in range(7):
            self.points.append(PlantThermoPoint(["Water"], [1]))

        # Point 8-9 - Re-injection section:
        #
        #   8 - Steam Cylinder Liquid Outlet
        #   9 - BHE Input

        self.points.append(PlantThermoPoint(["Water"], [1]))
        self.points.append(self.BHE_input)

        # Points 10 - 10 + 2(n - 1) - Inter-cooled Compression Section:
        # (2(n - 1) points - with n: number of stages)
        #
        #   10 + 2 * n       - Compression Outlet
        #   10 + 2 * n + 1   - Intercooler Outlet

        for i in range(self.__n_stages_comp - 1):
            self.points.append(PlantThermoPoint(["water"], [1]))
            self.points.append(PlantThermoPoint(["water"], [1]))

        # Points 10 + 2(n - 1) - end - Inter-cooler liquid handling Section:
        # (2(n - 1) + 2 points - with n: number of stages)
        # i_s_start is the section initial index: 10 + 2(n - 1) + 1
        #
        #   i_s_start               - Liquid intake for inter-cooling
        #   i_s_start + 2 * n       - n_th Pump Outlet
        #   i_s_start + 2 * n + 1   - n_th Intercooler Outlet
        #   i_end                   - cooling water discharge

        self.i_start_ic_water = len(self.points)
        self.points.append(PlantThermoPoint(["water"], [1]))

        for i in range(self.__n_stages_comp - 1):
            self.points.append(PlantThermoPoint(["water"], [1]))
            self.points.append(PlantThermoPoint(["water"], [1]))

        self.points.append(PlantThermoPoint(["water"], [1]))

        # properties of points 2 (saturated vapour), 4 (send condition), 5 (return condition), 8 (saturated liquid)
        # and i_s_start (saturated liquid) will remain fixed (only their flow rates will change)

        self.points[2].set_variable("T", self.T_sc)
        self.points[2].set_variable("x", 1)
        self.P_sc = self.points[2].get_variable("P")

        self.points[4].set_to_specific_super_heating(self.P_steam, self.IC_dt_SH)

        self.points[5].set_variable("P", self.P_steam)
        self.points[5].set_variable("x", 0)

        self.points[8].set_variable("T", self.T_sc)
        self.points[8].set_variable("x", 0)

        self.points[self.i_start_ic_water].set_variable("T", self.T_sc)
        self.points[self.i_start_ic_water].set_variable("x", 0)

    def __init_fixed_parameters(self):

        self.eta_comp = 0.8
        self.eta_pump = 0.8

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def BHE_fluid(self):

        return "water"

    def thermo_analysis(self):

        # self.BHE_HS=BHEHeatingSection(self.BHE_well, consider_pressure_losses=False)

        old_m_ratio_BHE = 0.
        for i in range(10):

            self.__thermo_analysis_step()

            err = abs(old_m_ratio_BHE - self.m_BHE) / self.m_BHE
            old_m_ratio_BHE = self.m_BHE
            if err < 10E-5:
                break

    def __thermo_analysis_step(self):

        # STEP - 1
        # Evaluate BHE input temperature and update well calculation

        P_BHE_in = self.BHE_input.get_variable("P")
        self.points[9].set_to_compression_result(P_BHE_in, self.eta_pump, self.points[8])
        self.BHE_well.update()

        # STEP - 2
        # Evaluate Steam cylinder input conditions

        h_in = self.points[0].get_variable("h")
        self.points[1].set_variable("h", h_in)
        self.points[1].set_variable("P", self.P_sc)

        # STEP - 3
        # Evaluate Inter-cooled Compression

        self.__evaluate_inter_cooled_compression()

        # STEP - 4
        # Evaluate return points

        dh_RH = self.points[3].get_variable("h") - self.points[4].get_variable("h")
        h_in = self.points[5].get_variable("h")

        self.points[6].set_variable("h", h_in + dh_RH)
        self.points[6].set_variable("P", self.P_steam)

        self.points[7].set_variable("h", self.points[6].get_variable("h"))
        self.points[7].set_variable("P", self.P_sc)

        self.__evaluate_flow_rates()
        self.__calculate_power()

    def __evaluate_inter_cooled_compression(self):

        # 1 - Initialization

        # list in which store power
        self.IC_pump_power_list = list()
        self.comp_power_list = list()
        self.ic_power_list = list()

        # overall coolant flow rate and output enthalpy
        self.points[self.i_start_ic_water].m_dot = 0.
        h_cool_out = 0.

        # 2 - Evaluation of compressor pressure ratio
        # Pressure Ratio is equal for each compressor:
        # beta_comp = beta_overall ^ (1 / n)
        beta_comp = np.power(self.P_steam / self.points[2].get_variable("P"), 1 / self.n_stages_comp)

        # 3 - Stages Iteration
        comp_in_i = 2
        for n in range(self.n_stages_comp - 1):
            # Steam Current Indices
            comp_out_i = 10 + 2 * n
            ic_out_i = comp_out_i + 1

            # Steam Point Calculation
            P_out = beta_comp * self.points[comp_in_i].get_variable("P")
            self.points[comp_out_i].set_to_compression_result(P_out, self.eta_comp, self.points[comp_in_i])
            self.points[ic_out_i].set_to_specific_super_heating(P_out, super_heating_value=self.IC_dt_SH)
            dt_ic = self.points[comp_out_i].get_variable("T") - self.points[ic_out_i].get_variable("T")

            # Cooling Water Current Indices
            pump_out_i = self.i_start_ic_water + 2 * n + 1
            ic_cool_out_i = pump_out_i + 1

            # Cooling Water Point Calculation
            T_out_coolant = self.points[self.i_start_ic_water].get_variable("T") + dt_ic
            P_out = self.__evaluate_cooling_water_pressure(ic_cool_out_i, T_out_coolant)
            self.points[pump_out_i].set_to_compression_result(P_out, self.eta_pump, self.points[self.i_start_ic_water])

            # Power Calculations
            dh_comp = self.points[comp_out_i].get_variable("h") - self.points[comp_in_i].get_variable("h")
            dh_pump = self.points[pump_out_i].get_variable("h") - self.points[self.i_start_ic_water].get_variable("h")

            dh_ic = self.points[comp_out_i].get_variable("h") - self.points[ic_out_i].get_variable("h")
            dh_ic_cool = self.points[ic_cool_out_i].get_variable("h") - self.points[pump_out_i].get_variable("h")
            m_cool = self.m_steam * dh_ic / dh_ic_cool
            h_cool_out += self.points[ic_cool_out_i].get_variable("h") * m_cool

            self.comp_power_list.append(dh_comp * self.m_steam)
            self.IC_pump_power_list.append(dh_pump * m_cool)
            self.ic_power_list.append(dh_ic * self.m_steam)

            # Set flow rates
            self.points[pump_out_i].m_dot = m_cool
            self.points[ic_cool_out_i].m_dot = m_cool
            self.points[self.i_start_ic_water].m_dot += m_cool

            # Update Steam Current Indices
            comp_in_i = ic_out_i

        # 4 - Last Stage Calculation
        self.points[3].set_to_compression_result(self.P_steam, self.eta_comp, self.points[comp_in_i])
        dh_comp = self.points[3].get_variable("h") - self.points[comp_in_i].get_variable("h")
        self.comp_power_list.append(dh_comp * self.m_steam)

        # 5 - Coolant Outlet Calculation
        m_cool_tot = self.points[self.i_start_ic_water].m_dot
        self.points[-1].set_variable("h", h_cool_out / m_cool_tot)
        self.points[-1].set_variable("P", self.P_sc)
        self.points[-1].m_dot = m_cool_tot

    def __evaluate_cooling_water_pressure(self, i, T_out_coolant):

        self.points[i].set_variable("T", T_out_coolant)
        self.points[i].set_variable("x", 0)

        P_sat = self.points[i].get_variable("P")
        self.points[i].set_variable("P", P_sat * 1.05)
        self.points[i].set_variable("T", T_out_coolant)

        return P_sat * 1.05

    def __evaluate_flow_rates(self):

        dh_steam = self.points[2].get_variable("h") - self.points[7].get_variable("h")
        dh_BHE = self.points[1].get_variable("h") - self.points[8].get_variable("h")
        dh_ic = self.points[-1].get_variable("h") - self.points[self.i_start_ic_water].get_variable("h")
        m_cool_tot = self.points[-1].m_dot

        self.m_BHE = self.m_steam * (dh_steam - dh_ic * m_cool_tot / self.m_steam) / dh_BHE
        self.m_dot_ratio = self.m_steam / self.m_BHE

        # set flow rates to points
        BHE_points = [0, 1, 8, 9]
        for i in range(len(self.points)):

            if i == self.i_start_ic_water:
                break

            if i in BHE_points:

                self.points[i].m_dot = self.m_BHE

            else:

                self.points[i].m_dot = self.m_steam

    def __calculate_power(self):

        m_dot_BHE = self.points[0].get_variable("m_dot")

        # steam produced
        dh_steam = self.points[4].get_variable("h") - self.points[5].get_variable("h")
        self.Q_steam = self.m_steam * dh_steam

        # BHE heat
        dh_BHE = self.points[0].get_variable("h") - self.points[9].get_variable("h")
        self.Q_BHE = m_dot_BHE * dh_BHE

        # inter_cooled_compression overall
        self.W_comp = sum(self.comp_power_list)
        self.W_IC = sum(self.ic_power_list)
        self.W_pump_IC = sum(self.IC_pump_power_list)

        # Re-heater
        dh_RH = self.points[3].get_variable("h") - self.points[4].get_variable("h")
        self.Q_RH = self.m_steam * dh_RH

        # pump
        dh_pump = self.points[9].get_variable("h") - self.points[8].get_variable("h")
        self.W_pump_BHE = m_dot_BHE * dh_pump

        self.W_net = self.W_comp + self.W_pump_BHE + self.W_pump_IC
        self.COP = self.Q_steam / self.W_net
        self.Q_ratio = self.Q_steam / self.Q_BHE

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update_economic_parameters(self):

        self.c_elet = 0.256 / 3600  # in €/kJ

    def calculate_PEC(self):

        self.comp_cost_list = list()
        self.pump_cost_list = list()
        self.he_area_list = list()
        self.ic_cost_list = list()

        self.RH_cost = self.__calculate_he_cost()
        self.BHE_pump_cost = TURTON_PEC_calculator(

            self.W_pump_BHE, 3.3892, 0.0536, 0.1538,
            B_1=1.89, B_2=1.35, F_M=2.2

        )

        for n in range(self.__n_stages_comp - 1):
            self.comp_cost_list.append(TURTON_PEC_calculator(

                self.comp_power_list[n], 2.2897, 1.3604, -0.1027,
                B_1=0, B_2=1, F_M=2.8

            ))

            self.pump_cost_list.append(TURTON_PEC_calculator(

                self.IC_pump_power_list[n], 3.3892, 0.0536, 0.1538,
                B_1=1.89, B_2=1.35, F_M=2.2

            ))

            self.ic_cost_list.append(self.__calculate_he_cost(n))

        self.comp_cost_tot = sum(self.comp_cost_list)
        self.ic_cost_tot = sum(self.ic_cost_list)
        self.pump_cost_tot = sum(self.pump_cost_list)

        return self.comp_cost_tot + self.ic_cost_tot + self.pump_cost_tot + self.RH_cost + self.BHE_pump_cost

    def __calculate_he_cost(self, n=None):

        # Heat Exchange Coefficient got from:
        # https://www.engineeringtoolbox.com/heat-transfer-coefficients-exchangers-d_450.html
        # U_ic defined from "Tubular, heating or cooling" and "Steam outside and liquid inside tubes"
        # U_RH higher as the incoming liquid is in phase changing state

        U_ic = 0.6  # in kW / (m^2 K)
        U_RH = 1  # in kW / (m^2 K)

        if n is None:

            # if n is none correlation of RH

            dt_a = abs(self.points[3].get_variable("T") - self.points[6].get_variable("T"))
            dt_b = abs(self.points[4].get_variable("T") - self.points[5].get_variable("T"))

            pressure = self.points[4].get_variable("P") * 10  # in bar
            C = [0.03881, -0.11272, 0.08183]  # from TURTON "pressure in shell and in tubes"

            dt = (dt_a - dt_b) / np.log(dt_a / dt_b)
            power = self.Q_RH
            U = U_RH

        else:

            comp_out_i = 10 + 2 * n
            ic_out_i = comp_out_i + 1
            pump_out_i = self.i_start_ic_water + 2 * n + 1

            pressure = self.points[pump_out_i].get_variable("P") * 10  # in bar
            C = [-0.00164, -0.00627, 0.0123]  # from TURTON "pressure in tubes only"

            dt = self.points[ic_out_i].get_variable("T") - self.points[self.i_start_ic_water].get_variable("T")
            power = self.ic_power_list[n]
            U = U_ic

        area = power / (U * dt)
        self.he_area_list.append(area)

        # TURTON Coefficient considering Fixed Tube Heat Exchanger:
        return TURTON_PEC_calculator(

            area, 4.3247, -0.3030, 0.1634, B_1=1.63, B_2=1.66, F_M=1.4,
            P=pressure, C_1=C[0], C_2=C[1], C_3=C[2]

        )

    def calculate_c_fuel(self):

        W_net = self.W_comp + self.W_pump_BHE
        eta_fuel = W_net / self.m_steam  # in kJ/kg

        return eta_fuel * self.c_elet

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------------   LCA METHODS   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def LCA_analysis(self):
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------------   OTHER METHODS   -------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def T_sc_perc(self):

        return self.__T_sc_perc

    @T_sc_perc.setter
    def T_sc_perc(self, T_sc_perc):

        self.__T_sc_perc = T_sc_perc
        self.__evaluate_T_sc()
        self.__init_points()

    @property
    def n_stages_comp(self):

        return self.__n_stages_comp

    @n_stages_comp.setter
    def n_stages_comp(self, n_stages_comp):

        self.__n_stages_comp = n_stages_comp
        self.__evaluate_T_sc()
        self.__init_points()


class WaterHeatPumpThermo(AbstractPythonHTHP):

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def optimization_bounds(self):

        if self.HTHP_fluid[0] == "Water":

            T_EVA_min = 2  # Supposed minimum water saturation temperature allowed (to avoid ice formation)

        else:

            T_EVA_min = -10  # Supposed minimum pentane saturation temperature allowed (to avoid ice formation)

        self.BHE_well.update()
        T_out = self.BHE_output.get_variable("T")
        range_max = T_out - T_EVA_min - self.DT_HE

        return Bounds([5], [range_max])

    @property
    def initial_guess(self):

        return [10]

    def set_optimization_param(self, optimization_param):

        self.range_eva = optimization_param[0]

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   INITIALIZATION METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def __init__(

            self,
            BHE_depth,
            T_rock,
            P_steam,

            m_steam=1.,
            T_in_BHE=30.,
            P_in_BHE=8.,
            T_ambient=15.,

            HTHP_fluid=None,
            HTHP_fluid_comp=None

    ):

        if HTHP_fluid is None:
            HTHP_fluid = ["Pentane"]

        if HTHP_fluid_comp is None:
            HTHP_fluid_comp = np.ones(len(HTHP_fluid))
            HTHP_fluid_comp = HTHP_fluid_comp / np.sum(HTHP_fluid_comp)

        self.HTHP_fluid = HTHP_fluid
        self.HTHP_fluid_comp = HTHP_fluid_comp

        super().__init__(BHE_depth, T_rock, P_steam, m_steam, T_in_BHE, P_in_BHE, T_ambient)

    def initialize_other_parameter(self):

        self.DT_HE = 15  # Pinch point in the EVA and COND heat exchangers [°C]
        self.range_eva = 10  # condenser water temperature range [°C] (To be optimized with economic considerations)
        self.eta_RH = 0.90  # Regenerator Efficiency
        self.eta_comp = 0.85  # Compressor Efficiency

        self.__init_points()

    def __init_points(self):

        self.points = list()

        # Point 0-1 - BHE Output/Input
        #
        #   0 - BHE Output
        #   1 - BHE Input

        self.points.append(self.BHE_output)
        self.points.append(self.BHE_input)

        # Points 2-7 - steam distribution section:
        #
        #   2 - EVA Outlet
        #   3 - Compressor Inlet
        #   4 - Compressor Outlet
        #   5 - COND Outlet
        #   6 - RH Outlet (hot side)
        #   7 - EVA inlet

        for i in range(6):
            self.points.append(PlantThermoPoint(self.HTHP_fluid, self.HTHP_fluid_comp))

        # Point 8-9 - Steam Production Section:
        #
        #   8 - Steam Cylinder Liquid Phase
        #   9 - Steam Cylinder Vapour Phase

        for i in range(2):
            self.points.append(PlantThermoPoint(["Water"], [1]))

        # properties of points 5 (fluid saturated liquid), 8 (saturated water), 9 (saturated steam)
        # will remain fixed (only their flow rates will change)

        self.points[8].set_variable("P", self.P_steam)
        self.points[8].set_variable("x", 0)

        self.points[9].set_variable("P", self.P_steam)
        self.points[9].set_variable("x", 1)

        self.points[5].set_variable("T", self.points[9].get_variable("T") + self.DT_HE)
        self.points[5].set_variable("x", 0)

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def BHE_fluid(self):
        return "water"

    def thermo_analysis(self):

        # STEP - 1
        # Evaluate BHE input temperature and update well calculation

        self.BHE_input.set_variable("T", self.BHE_output.get_variable("T") - self.range_eva)
        self.BHE_well.update()

        # STEP - 2
        # Evaluate COND outlet conditions

        self.points[2].set_variable("T", self.BHE_input.get_variable("T") - self.DT_HE)
        self.points[2].set_variable("x", 1)

        # STEP - 3
        # Evaluate RH conditions

        self.__evaluate_RH()

        # STEP - 4
        # Evaluate Compression and Expansion
        self.points[4].set_to_compression_result(

            P_out=self.points[5].get_variable("P"),
            input_state=self.points[3],
            eta=self.eta_comp

        )

        self.points[7].set_variable("P", self.points[2].get_variable("P"))
        self.points[7].set_variable("h", self.points[6].get_variable("h"))

        self.__evaluate_flow_rates()
        self.__calculate_power()

        # print()
        # print("{} {} -> {} {}".format(self.range_eva, self.eta_RH, self.COP, self.m_BHE))

    def __evaluate_RH(self):

        # STEP - 0
        # Check regeneration feasibility

        if self.points[2].get_variable("T") > self.points[5].get_variable("T"):
            # if temperature in point 2 is higher than temperature in point 5 regeneration makes no sense,
            # hence no heat exchange is expected and temperatures remains the same

            self.points[3].set_variable("P", self.points[2].get_variable("P"))
            self.points[6].set_variable("P", self.points[5].get_variable("P"))

            self.points[3].set_variable("T", self.points[2].get_variable("T"))
            self.points[6].set_variable("T", self.points[5].get_variable("T"))

            return

        # STEP - 1
        # Identify max heat transfer

        self.points[3].set_variable("P", self.points[2].get_variable("P"))
        self.points[6].set_variable("P", self.points[5].get_variable("P"))

        self.points[3].set_variable("T", self.points[5].get_variable("T"))
        self.points[6].set_variable("T", self.points[2].get_variable("T"))

        HE_max_cold = self.points[3].evaluate_variable_variation(self.points[2], "h")
        HE_max_hot = self.points[5].evaluate_variable_variation(self.points[6], "h")
        HE_max = min(abs(HE_max_cold), abs(HE_max_hot))
        HE = HE_max * self.eta_RH

        # STEP - 2
        # Evaluate Results

        self.points[3].set_variable("P", self.points[2].get_variable("P"))
        self.points[3].set_variable("h", self.points[2].get_variable("h") + HE)

        self.points[6].set_variable("P", self.points[5].get_variable("P"))
        self.points[6].set_variable("h", self.points[5].get_variable("h") - HE)

    def __evaluate_flow_rates(self):

        dh_steam = self.points[9].evaluate_variable_variation(self.points[8], "h")
        dh_COND = self.points[4].evaluate_variable_variation(self.points[5], "h")
        dh_EVA = self.points[2].evaluate_variable_variation(self.points[7], "h")
        dh_BHE_brine = self.points[0].evaluate_variable_variation(self.points[1], "h")

        self.m_ORC = dh_steam / dh_COND * self.m_steam
        self.m_BHE = self.m_ORC * dh_EVA / dh_BHE_brine
        self.m_dot_ratio = self.m_steam / self.m_BHE
        self.dh_ratio = dh_steam/dh_BHE_brine

        for i in range(2):
            self.points[i].m_dot = self.m_BHE

        for i in range(6):
            self.points[i + 2].m_dot = self.m_ORC

        for i in range(2):
            self.points[i + 6 + 2].m_dot = self.m_steam

    def __calculate_power(self):

        dh_BHE = self.points[0].evaluate_variable_variation(self.points[1], "h")
        dh_RH = self.points[3].evaluate_variable_variation(self.points[2], "h")
        dh_COND = self.points[4].evaluate_variable_variation(self.points[5], "h")
        dh_COMP = self.points[4].evaluate_variable_variation(self.points[3], "h")

        self.Q_steam = self.m_ORC * dh_COND
        self.Q_BHE = self.m_BHE * dh_BHE
        self.Q_RH = self.m_ORC * dh_RH
        self.W_COMP = self.m_ORC * dh_COMP
        self.Q_COND = self.Q_steam
        self.Q_EVA = self.Q_BHE

        self.COP = self.Q_steam / self.W_COMP
        self.Q_ratio = self.Q_steam / self.Q_BHE

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def calculate_PEC(self):
        return 0.

    def calculate_c_fuel(self):
        return 0.

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------------   LCA METHODS   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def LCA_analysis(self):
        pass
