from main_code.power_plants.HTHP.abstract_hthp import AbstractPythonHTHP
from main_code.support import PlantThermoPoint
from scipy.optimize import Bounds
import numpy as np


class CO2HeatPumpThermo(AbstractPythonHTHP):

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def optimization_bounds(self):

        return Bounds([0, 700], [1, 900])
        # return None

    @property
    def initial_guess(self):

        return [0.5, 800]

    def set_optimization_param(self, optimization_param):

        self.T_SG_perc = optimization_param[0]
        self.rho_in = optimization_param[1]

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   INITIALIZATION METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def initialize_other_parameter(self):

        self.T_SG_perc = 0.5  # Value to be optimized that control the temperature of the sCO2 (Limited Between 0-1)
        self.dT_SG_pinch = 8  # Temperature difference between water and sCO2 in the SG heat exchanger [°C]
        self.P_CO2_max = 100  # Maximum CO2 pressure allowed in the system [MPa] (overestimated to extend the optimization range)
        self.rho_in = 800     # BHE input density [kg/m3]

        self.__init_fixed_parameters()
        self.__init_points()
        self.__init_economic_results()

    def __init_fixed_parameters(self):

        self.eta_comp = 0.8
        self.eta_turb = 0.75
        self.eta_pump = 0.8

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

        self.points[5].set_variable("P", self.P_steam)
        self.points[5].set_variable("x", 0)

        self.points[9].set_variable("P", self.P_steam)
        self.points[9].set_variable("x", 1)

    def __init_economic_results(self):

        self.c_HTHE = 0.
        self.c_GC = 0.
        self.c_turb = 0.
        self.c_comp = 0.
        self.c_motor = 0.

        return self.PEC

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def thermo_analysis(self):

        self.BHE_input.set_variable("T", self.BHE_input.get_variable("T"))
        self.BHE_input.set_variable("rho", self.rho_in)

        self.BHE_well.update()

        # STEP - 1
        # Calculate T_CO_2_max - Temperature reached by the CO2 if the maximum CO2 pressure is reached
        # (Point 1 used as a abstract_classes point)

        self.points[1].set_to_compression_result(self.P_CO2_max, self.eta_comp, input_state=self.points[0])
        T_CO_2_max = self.points[1].get_variable("T")

        # T_SG_HE_out evaluation - the value should remain between the minimum between T_steam and T_BHE_out (minimum)
        # and T_CO_2_max - DT_HE (maximum) - T_SG_perc identify the position in between the two

        T_SG_HE_out_min = max(self.points[5].get_variable("T"), self.points[0].get_variable("T"))
        T_SG_HE_out_max = T_CO_2_max - self.dT_SG_pinch
        T_SG_HE_out = (T_SG_HE_out_max - T_SG_HE_out_min) * self.T_SG_perc + T_SG_HE_out_min

        # STEP - 2
        # P_SG_HE evaluation (5% more than the saturation pressure for the selected temperature or, in case, more than
        # the critical pressure, in order to assure liquid condition) - Point 6 used as a abstract_classes point

        if T_SG_HE_out < self.points[5].RPHandler.TC:

            self.points[6].set_variable("T", T_SG_HE_out)
            self.points[6].set_variable("x", 0)
            P_SG_HE = self.points[6].get_variable("P") * 1.05

        else:

            P_SG_HE = self.points[5].RPHandler.PC * 1.05

        # STEP - 3
        # Steam Point Calculation

        # Point 6 - Circulation Pump Output
        self.points[6].set_to_compression_result(P_SG_HE, self.eta_pump, input_state=self.points[5])

        # Point 7 - Water Main HE Output
        self.points[7].set_variable("T", T_SG_HE_out)
        self.points[7].set_variable("P", P_SG_HE)

        # Point 8 - Lamination Valve Output
        self.points[8].set_variable("h", self.points[7].get_variable("h"))
        self.points[8].set_variable("P", self.P_steam)

        # STEP - 4
        # point 1 iteration

        T_CO2_max = T_SG_HE_out + self.dT_SG_pinch
        self.__iterate_P_max(T_CO2_max)

        # STEP - 5
        # CO2 Point Calculation

        # Point 2 - CO2 Main HE Output
        self.points[2].set_variable("T", self.points[6].get_variable("T") + self.dT_SG_pinch)
        self.points[2].set_variable("P", self.points[1].get_variable("P"))

        # Point 3 - CO2 Expander Output
        self.points[3].set_to_expansion_result(self.points[4].get_variable("P"), self.eta_turb, self.points[2])

        self.__evaluate_flow_rates()
        self.__calculate_power()

    def __iterate_P_max(self, T_CO2_HE):

        # iteration uses bisection method for speed, P_up initial is P_CO2_max P_down initial is P_BHE_out

        toll = 0.001
        error = 10 * toll
        self.iter = 0
        P_limit = [self.points[0].get_variable("P"), self.P_CO2_max]

        while error > toll:

            P_mean = sum(P_limit) / 2
            self.points[1].set_to_compression_result(P_mean, self.eta_comp, self.points[0])

            T_res = self.points[1].get_variable("T")

            if T_res > T_CO2_HE:

                P_limit[1] = P_mean

            else:

                P_limit[0] = P_mean

            error = abs(T_res - T_CO2_HE) / T_CO2_HE
            self.iter += 1
            # print(error)

            if self.iter > 200:
                break

    def __evaluate_flow_rates(self):

        steam_water_flow_ratio = self.points[8].get_variable("x")

        dh_CO2 = self.points[1].get_variable("h") - self.points[2].get_variable("h")
        dh_water = self.points[7].get_variable("h") - self.points[6].get_variable("h")

        self.steam_CO2_flow_ratio = dh_CO2 / dh_water * steam_water_flow_ratio
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

    # <------------------------------------------------------------------------->

    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update_economic_parameters(self):

        self.c_elet = 0.256 / 3600  # in €/kJ

    def calculate_PEC(self):

        self.PEC = 0.
        m_dot_co2 = self.points[0].get_variable("m_dot")

        # High Temperature Heat Exchanger Cost
        UA_HTHE = self.__calculate_UA_HTHE()
        self.c_HTHE = 49.45 * np.power(UA_HTHE, 0.7544)
        self.PEC += self.c_HTHE

        # Gas Cooler Cost
        UA_GC = self.calculate_UA_GC()
        self.c_GC = 32.88 * np.power(UA_GC, 0.75)
        self.PEC += self.c_GC

        # Turbine Cost
        self.c_turb = 406200 * np.power(self.W_turb / 1000, 0.8)  # in MW
        self.PEC += self.c_turb

        # Compressor Cost
        self.c_comp = 1230000 * np.power(self.W_comp / 1000, 0.3992)  # in MW
        self.PEC += self.c_comp

        # Motor (or generator) Cost
        self.c_motor = 211400 * np.power(abs(self.W_comp - self.W_turb) / 1000, 0.6227)  # in MW
        self.PEC += self.c_motor

        return self.PEC

    def calculate_c_fuel(self):

        W_net = self.W_comp + self.W_pump - self.W_turb
        eta_fuel = W_net / self.m_steam  # in kJ/kg

        return eta_fuel * self.c_elet

    def __calculate_UA_HTHE(self):

        LMTD = self.dT_SG_pinch
        return self.Q_HTHE * 1000 / LMTD  # in W/K

    def calculate_UA_GC(self):

        DT_A = self.points[3].get_variable("T") - self.ambient.get_variable("T")
        DT_B = self.points[4].get_variable("T") - self.ambient.get_variable("T")

        LMTD = (DT_A - DT_B) / np.log(DT_A / DT_B)
        return self.Q_GC * 1000 / LMTD  # in W/K

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------------   LCA METHODS   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def LCA_analysis(self):
        pass


class CO2HeatPumpThermoRegeneration(AbstractPythonHTHP):

    # TODO

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def optimization_bounds(self):

        return Bounds([0, 700], [1, 900])
        # return None

    @property
    def initial_guess(self):

        return [0.5, 800]

    def set_optimization_param(self, optimization_param):

        self.T_SG_perc = optimization_param[0]
        self.rho_in = optimization_param[1]

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   INITIALIZATION METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def initialize_other_parameter(self):

        self.T_SG_perc = 0.5  # Value to be optimized that control the temperature of the sCO2 (Limited Between 0-1)
        self.dT_SG_pinch = 8  # Temperature difference between water and sCO2 in the SG heat exchanger [°C]
        self.P_CO2_max = 100  # Maximum CO2 pressure allowed in the system [MPa] (overestimated to extend the optimization range)
        self.rho_in = 800     # BHE input density [kg/m3]

        self.__init_fixed_parameters()
        self.__init_points()
        self.__init_economic_results()

    def __init_fixed_parameters(self):

        self.eta_comp = 0.8
        self.eta_turb = 0.75
        self.eta_pump = 0.8

        self.eta_RH = 0.90  # Regenerator Efficiency

    def __init_points(self):

        self.points = list()

        # Point 0 - BHE Output

        self.points.append(self.BHE_output)

        # Points 1-5 - CO2 BHE 2022-10-04 - HTHP internal points:
        #
        #   1 - CO2 Cold RH Output
        #   2 - CO2 Compressor Output
        #   3 - CO2 Main HE Output
        #   4 - CO2 Hot RH Output
        #   5 - CO2 Expander Output

        for i in range(5):
            self.points.append(PlantThermoPoint(["CarbonDioxide"], [1]))

        # Point 6 - BHE Input
        self.points.append(self.BHE_input)

        # Points 7-11 - Steam Drum internal points:
        #
        #   07 - Steam Drum Liquid Output
        #   08 - Circulation Pump Output
        #   09 - Water Main HE Output
        #   10 - Lamination Valve Output
        #   11 - Steam Drum Vapour Output

        for i in range(5):
            self.points.append(PlantThermoPoint(["water"], [1]))

        # point 7 (saturated liquid) and 11 (saturated vapour) properties will remain fixed (only their flow rates
        # will change)

        self.points[7].set_variable("P", self.P_steam)
        self.points[7].set_variable("x", 0)

        self.points[11].set_variable("P", self.P_steam)
        self.points[11].set_variable("x", 1)

    def __init_economic_results(self):

        self.c_HTHE = 0.
        self.c_GC = 0.
        self.c_turb = 0.
        self.c_comp = 0.
        self.c_motor = 0.

        return self.PEC

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def thermo_analysis(self):

        self.BHE_input.set_variable("T", self.BHE_input.get_variable("T"))
        self.BHE_input.set_variable("rho", self.rho_in)
        self.BHE_well.update()

        # STEP - 1
        # Evaluate RH

        self.__evaluate_RH()

        # STEP - 2
        # Calculate T_CO_2_max - Temperature reached by the CO2 if the maximum CO2 pressure is reached
        # (Point 1 used as a abstract_classes point)

        self.points[1].set_to_compression_result(self.P_CO2_max, self.eta_comp, input_state=self.points[0])
        T_CO_2_max = self.points[1].get_variable("T")

        # T_SG_HE_out evaluation - the value should remain between:
        #
        #   (minimum) -> the minimum between T_steam and T_RH_out
        #   (maximum) -> T_CO_2_max - DT_HE
        #
        #   T_SG_perc identify the position in between the two extremes

        T_SG_HE_out_min = max(self.points[5].get_variable("T"), self.points[0].get_variable("T"))
        T_SG_HE_out_max = T_CO_2_max - self.dT_SG_pinch
        T_SG_HE_out = (T_SG_HE_out_max - T_SG_HE_out_min) * self.T_SG_perc + T_SG_HE_out_min

        # STEP - 2
        # P_SG_HE evaluation (5% more than the saturation pressure for the selected temperature or, in case, more than
        # the critical pressure, in order to assure liquid condition) - Point 6 used as a abstract_classes point

        if T_SG_HE_out < self.points[5].RPHandler.TC:

            self.points[6].set_variable("T", T_SG_HE_out)
            self.points[6].set_variable("x", 0)
            P_SG_HE = self.points[6].get_variable("P") * 1.05

        else:

            P_SG_HE = self.points[5].RPHandler.PC * 1.05

        # STEP - 3
        # Steam Point Calculation

        # Point 6 - Circulation Pump Output
        self.points[6].set_to_compression_result(P_SG_HE, self.eta_pump, input_state=self.points[5])

        # Point 7 - Water Main HE Output
        self.points[7].set_variable("T", T_SG_HE_out)
        self.points[7].set_variable("P", P_SG_HE)

        # Point 8 - Lamination Valve Output
        self.points[8].set_variable("h", self.points[7].get_variable("h"))
        self.points[8].set_variable("P", self.P_steam)

        # STEP - 4
        # point 1 iteration

        T_CO2_max = T_SG_HE_out + self.dT_SG_pinch
        self.__iterate_P_max(T_CO2_max)

        # STEP - 5
        # CO2 Point Calculation

        # Point 2 - CO2 Main HE Output
        self.points[2].set_variable("T", self.points[6].get_variable("T") + self.dT_SG_pinch)
        self.points[2].set_variable("P", self.points[1].get_variable("P"))

        # Point 3 - CO2 Expander Output
        self.points[3].set_to_expansion_result(self.points[4].get_variable("P"), self.eta_turb, self.points[2])

        self.__evaluate_flow_rates()
        self.__calculate_power()

    def __evaluate_RH(self):

        # STEP - 0
        # Identify hot side input temperature

        if self.points[2].get_variable("T") > self.points[5].get_variable("T"):

            # if temperature in point 0 is higher than temperature in point 3 regeneration makes no sense,
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

        HE_max_cold = self.points[3].dvar(self.points[2], "h")
        HE_max_hot = self.points[5].dvar(self.points[6], "h")
        HE_max = min(abs(HE_max_cold), abs(HE_max_hot))
        HE = HE_max * self.eta_RH

        # STEP - 2
        # Evaluate Results

        self.points[3].set_variable("P", self.points[2].get_variable("P"))
        self.points[3].set_variable("h", self.points[2].get_variable("h") + HE)

        self.points[6].set_variable("P", self.points[5].get_variable("P"))
        self.points[6].set_variable("h", self.points[5].get_variable("h") - HE)

    def __iterate_P_max(self, T_CO2_HE):

        # iteration uses bisection method for speed, P_up initial is P_CO2_max P_down initial is P_BHE_out

        toll = 0.001
        error = 10 * toll
        self.iter = 0
        P_limit = [self.points[0].get_variable("P"), self.P_CO2_max]

        while error > toll:

            P_mean = sum(P_limit) / 2
            self.points[1].set_to_compression_result(P_mean, self.eta_comp, self.points[0])

            T_res = self.points[1].get_variable("T")

            if T_res > T_CO2_HE:

                P_limit[1] = P_mean

            else:

                P_limit[0] = P_mean

            error = abs(T_res - T_CO2_HE) / T_CO2_HE
            self.iter += 1
            # print(error)

            if self.iter > 200:
                break

    def __evaluate_flow_rates(self):

        steam_water_flow_ratio = self.points[8].get_variable("x")

        dh_CO2 = self.points[1].get_variable("h") - self.points[2].get_variable("h")
        dh_water = self.points[7].get_variable("h") - self.points[6].get_variable("h")

        self.steam_CO2_flow_ratio = dh_CO2 / dh_water * steam_water_flow_ratio
        self.m_BHE = self.m_steam / self.steam_CO2_flow_ratio
        self.m_dot_ratio = self.m_steam / self.m_BHE

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

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update_economic_parameters(self):

        self.c_elet = 0.256 / 3600  # in €/kJ

    def calculate_PEC(self):

        self.PEC = 0.
        m_dot_co2 = self.points[0].get_variable("m_dot")

        # High Temperature Heat Exchanger Cost
        UA_HTHE = self.__calculate_UA_HTHE()
        self.c_HTHE = 49.45 * np.power(UA_HTHE, 0.7544)
        self.PEC += self.c_HTHE

        # Gas Cooler Cost
        UA_GC = self.__calculate_UA_GC()
        self.c_GC = 32.88 * np.power(UA_GC, 0.75)
        self.PEC += self.c_GC

        # Turbine Cost
        self.c_turb = 406200 * np.power(self.W_turb / 1000, 0.8)  # in MW
        self.PEC += self.c_turb

        # Compressor Cost
        self.c_comp = 1230000 * np.power(self.W_comp / 1000, 0.3992)  # in MW
        self.PEC += self.c_comp

        # Motor (or generator) Cost
        self.c_motor = 211400 * np.power(abs(self.W_comp - self.W_turb) / 1000, 0.6227)  # in MW
        self.PEC += self.c_motor

        return self.PEC

    def calculate_c_fuel(self):

        W_net = self.W_comp - self.W_turb
        eta_fuel = W_net / self.m_steam  # in kJ/kg

        return eta_fuel * self.c_elet

    def __calculate_UA_HTHE(self):

        LMTD = self.dT_SG_pinch
        return self.Q_HTHE * 1000 / LMTD  # in W/K

    def __calculate_UA_GC(self):

        DT_A = self.points[3].get_variable("T") - self.ambient.get_variable("T")
        DT_B = self.points[4].get_variable("T") - self.ambient.get_variable("T")

        LMTD = (DT_A - DT_B) / np.log(DT_A / DT_B)
        return self.Q_GC * 1000 / LMTD  # in W/K

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------------   LCA METHODS   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def LCA_analysis(self):
        pass


class StandaloneCO2HeatPump(CO2HeatPumpThermo):

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   INITIALIZATION METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def __append_other_points(self):

        # Point 0-9 - Standard sCO2 2022-10-04 - HTHP
        # Point 10 - Turbine Outlet
        # Point 11 - Mixer before GC

        self.points.append(PlantThermoPoint(["CarbonDioxide"], [1]))
        self.points.append(PlantThermoPoint(["CarbonDioxide"], [1]))

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def thermo_analysis(self):

        self.__append_other_points()

        # STEP - 1
        # Calculate standard 2022-10-04 - HTHP
        super(StandaloneCO2HeatPump, self).thermo_analysis()

        # STEP - 2
        # Evaluate Turbine output state
        self.points[10].set_to_expansion_result(self.points[4].get_variable("P"), self.eta_turb, self.points[0])

        # STEP - 3
        # Update flow rates
        self.__update_flow_rates()

        # STEP - 4
        # Evaluate GC input state
        h_GC_in = (self.m_dot_HTHP * self.points[3].get_variable("h") + self.m_dot_turb_LP * self.points[10].get_variable("h")) / self.m_BHE
        self.points[11].set_variable("P", self.points[4].get_variable("P"))
        self.points[11].set_variable("h", h_GC_in)

        # STEP - 5
        # Update Power calculation
        self.__update_power()

    def __update_flow_rates(self):

        self.W_net_HTHP = self.W_comp + self.W_pump - self.W_turb

        if self.W_net_HTHP > 0:

            dh_turb_LP = self.points[0].get_variable("h") - self.points[10].get_variable("h")
            self.m_dot_turb_LP = self.W_net_HTHP / dh_turb_LP

        else:

            self.m_dot_turb_LP = 0.

        self.m_dot_ratio_HTHP = self.m_dot_ratio
        self.m_dot_HTHP = self.m_BHE

        self.m_BHE = self.m_dot_HTHP + self.m_dot_turb_LP
        self.m_dot_ratio_real = self.m_steam / self.m_BHE
        self.m_dot_ratio = 2.45 * self.m_steam / self.m_BHE

        # set flow rates to points

        # BHE points
        self.points[0].m_dot = self.m_BHE
        self.points[4].m_dot = self.m_BHE
        self.points[11].m_dot = self.m_BHE

        # LP Turbine points
        self.points[10].m_dot = self.m_dot_turb_LP

    def __update_power(self):

        # Gas Cooler
        dh_GC = self.points[11].get_variable("h") - self.points[4].get_variable("h")
        self.Q_GC = dh_GC * self.m_BHE

        # BHE
        dh_BHE = self.BHE_output.get_variable("h") - self.BHE_input.get_variable("h")
        self.Q_BHE = self.m_BHE * dh_BHE

        # TURB_LP
        dh_turb_LP = self.points[0].get_variable("h") - self.points[10].get_variable("h")
        self.W_turb_LP = self.m_dot_turb_LP * dh_turb_LP

        self.Q_ratio = self.Q_HTHE / self.Q_BHE
        self.COP = 0

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def calculate_PEC(self):

        super(StandaloneCO2HeatPump, self).calculate_PEC()

        # Turbine LP Cost
        self.c_turb_LP = 406200 * np.power(self.W_turb_LP / 1000, 0.8)  # in MW
        self.PEC += self.c_turb_LP

        return self.PEC

    def calculate_c_fuel(self):

        return 0.

    def calculate_UA_GC(self):

        DT_A = self.points[11].get_variable("T") - self.ambient.get_variable("T")
        DT_B = self.points[4].get_variable("T") - self.ambient.get_variable("T")

        LMTD = (DT_A - DT_B) / np.log(DT_A / DT_B)
        return self.Q_GC * 1000 / LMTD  # in W/K

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------------   LCA METHODS   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def LCA_analysis(self):
        pass
