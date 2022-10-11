from main_code.power_plants.abstract_surface_plant import AbstractSurfacePlant
from main_code.simplified_BHE.simplified_BHE import SimplifiedBHE
from main_code.support import PlantThermoPoint
from abc import ABC, abstractmethod
from scipy.optimize import Bounds
import scipy.optimize
import numpy as np


class AbstractPythonHTHP(AbstractSurfacePlant, ABC):

    def __init__(

            self,

            BHE_depth,
            T_rock,

            P_steam,
            m_steam=1,

            T_in_BHE=30,
            P_in_BHE=8,
            T_ambient=15,

            HTHP_fluid=None,
            HTHP_fluid_comp=None

    ):

        super().__init__()

        self.P_steam = P_steam
        self.m_steam = m_steam

        # main results
        self.COP = 0.
        self.Q_ratio = 0.
        self.m_dot_ratio = 0.

        # optimization parameters
        self.optimize_economic = False
        self.omega = 1.
        self.optimization_result = list()

        self.__init_BHE_well(BHE_depth, T_rock, T_in_BHE, P_in_BHE, T_ambient)
        self.__init_economic_parameters()
        self.initialize_other_parameter()

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   INITIALIZATION METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def BHE_fluid(self):

        """

            This property can be overwritten in sub-classes and used to change the working point to be used in the BHE
            the default is "CarbonDioxide" and will be returned if the property is not overwritten

        """
        return "CarbonDioxide"

    def __init_BHE_well(self, BHE_depth, T_rock, T_in_BHE, P_in_BHE, T_ambient):

        BHE_fluid = self.BHE_fluid
        input_point = PlantThermoPoint([BHE_fluid], [1])
        ambient_point = PlantThermoPoint([BHE_fluid], [1])

        input_point.set_variable("T", T_in_BHE)
        input_point.set_variable("P", P_in_BHE)

        ambient_point.set_variable("P", 0.1)
        ambient_point.set_variable("T", T_ambient)

        self.BHE_well = SimplifiedBHE(input_point, dz_well=BHE_depth, T_rocks=T_rock)
        self.append_ambient_condition(ambient_point)

        self.BHE_well.update(ambient_point)
        self.append_BHE_well_points(

            BHE_input=self.BHE_well.points[0],
            BHE_output=self.BHE_well.points[-1]

        )

    def __init_economic_parameters(self):

        self.c_steam = 0.
        self.PEC = 0.
        self.c_well = 0.
        self.c_fuel = 0.

        self.n_life = 20
        self.i_rate = 0.04
        self.t_op = 8000 * 3600  # in seconds
        self.gamma = 0.055
        self.beta = 2.40

    def initialize_other_parameter(self):

        """

            This function can be overwritten in sub-classes and used to initialize support parameter without overwriting
            the __init__ function (it has a lot of parameters and the code will look cleaner)

        """
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <--------------------   OPTIMIZATION METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def optimize_thermo(self, omega=1):

        self.omega = omega
        bounds = self.optimization_bounds

        if bounds is None:

            res = scipy.optimize.minimize(

                self.optimization_function,
                self.initial_guess

            )

        else:

            res = scipy.optimize.minimize(

                self.optimization_function,
                self.initial_guess,
                bounds=bounds

            )

        self.optimization_result = res.x
        self.set_optimization_param(res.x)
        self.calculate(calculate_thermo_only=True)

    def optimization_function(self, optimization_param):

        self.set_optimization_param(optimization_param)
        self.calculate(calculate_thermo_only=(not self.optimize_economic))
        return 1 / (self.m_dot_ratio * 20) + self.omega * 1 / (self.COP / 3)

    @property
    def overall_initial_guess(self):

        overall_initial_guess = list()
        overall_initial_guess.extend(self.initial_guess)
        overall_initial_guess.extend(self.BHE_well.initial_guess)
        return overall_initial_guess

    @property
    def overall_optimization_bounds(self):

        overall_initial_bound = self.optimization_bounds
        lb = np.append(overall_initial_bound.lb, self.BHE_well.optimization_bounds.lb)
        ub = np.append(overall_initial_bound.ub, self.BHE_well.optimization_bounds.ub)
        return Bounds(lb, ub)

    def set_overall_optimization_param(self, optimization_param):

        n_param_self = len(self.initial_guess)
        self.set_optimization_param(optimization_param[:n_param_self])
        self.BHE_well.set_optimization_param(optimization_param[n_param_self:])

    @property
    @abstractmethod
    def initial_guess(self):
        """
        This function must be overwritten in sub-classes and must be used to return the first guess for the optimization
        process
        """
        pass

    @property
    @abstractmethod
    def optimization_bounds(self) -> Bounds:
        """
        This function must be overwritten in sub-classes and must be used to return the optimization bonds for the plant
        """

    @abstractmethod
    def set_optimization_param(self, optimization_param):
        """
        This function must be overwritten in sub-classes and must be used to set the optimization parameter in
        the sub_class
        """
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-----------------   ECONOMIC ANALYSIS METHODS   ------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def economic_analysis(self):

        self.update_economic_parameters()
        self.BHE_well.update_economic()
        alpha = self.__calculate_alpha()

        self.PEC = self.calculate_PEC()
        self.c_aux = self.calculate_auxiliary_costs()
        self.c_well = self.BHE_well.c_well
        self.c_fuel = self.calculate_c_fuel()

        self.c_tot = self.beta * self.PEC + self.c_well + self.c_aux
        self.c_steam = (1 + self.gamma * alpha) / (self.m_steam * self.t_op) * self.c_tot + self.c_fuel

    def __calculate_alpha(self):

        alpha = 1 / self.i_rate * (1 - np.power((1 + self.i_rate), -self.n_life))
        return alpha

    def update_economic_parameters(self):

        """

            This function can be overwritten in sub-classes and can overwrite the value of:

                - self.n_life = 20, the expected operational lifetime of the plant in years
                - self.i_rate = 0.04, the expected yearly interest rate
                - self.t_op = 8000 * 3600, the anticipated operational time of the plant in a year (s/year)
                - self.gamma = 0.055, Ratio between the investment cost and the annual O&M cost
                - self.beta = 2.4, Ratio between the investment cost and the primary equipment cost

            if the function is not overwritten, default values shown above will be considered for the analysis

        """
        pass

    @abstractmethod
    def calculate_PEC(self):

        """
            This function must be overwritten in sub-classes and must return the PEC of the surface plant in €
        """
        pass

    @abstractmethod
    def calculate_c_fuel(self):
        """
        This function must be overwritten in sub-classes and must return the sum of the cost of alla the fuels needed
        for the plant to work in € (electrical power etc.)
        """
        pass

    def calculate_auxiliary_costs(self):

        """

            This method can be overwritten in sub-classes and can be used to add support auxiliary cost to the plant
            For example the cost the solar fields installed to feed the system

        """
        return 0.