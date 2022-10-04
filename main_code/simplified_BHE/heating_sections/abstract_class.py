from abc import ABC, abstractmethod
import warnings


class AbstractHeatingSection(ABC):

    def __init__(self, main_BHE):

        self.main_BHE = main_BHE
        self.main_BHE.heating_section = self
        self.condition_changed = True

        self.l_HS = 0.
        self.dP_HS = 0.
        self.CAPEX_well = 0.

    def __setattr__(self, attr, value):

        """

            The set attribute defined in this way means that once a private attribute (i.e. a geometric parameter) is
            modified, "condition_changed" become False, signifying that the thermo model must be recalculated.

        """

        object.__setattr__(self, attr, value)

        if "condition_changed" not in attr and type(self).__name__ in attr:

            self.condition_changed = True

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------   PROPERTIES   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def input_point(self):

        if len(self.main_BHE.points) > 1:

            return self.main_BHE.points[1]

        else:

            return None

    @property
    def output_point(self):

        if len(self.main_BHE.points) > 2:

            return self.main_BHE.points[2]

        else:

            return None

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   THERMO ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def update(self):

        if self.input_point is not None:

            if self.output_point is None:

                self.main_BHE.points.append(self.input_point.duplicate())

            self.update_HS()
            self.condition_changed = False

    @abstractmethod
    def update_HS(self):

        """

        This method must be updated in sub-classes, is called during the calculation of the simplified well and must
        update the value of "output_point" given the condition of "input_thermo_point" and the geometry of the well

        """

        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   ECONOMIC ANALYSIS METHODS   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def evaluate_cost(self):

        """

            Return the expected investment cost for well drilling [€]
            the most important parameter to consider are:

                 - dT_rocks [°C], Temperature difference between fluid output and undisturbed rocks

                 - n_wells [-], number of wells to be drilled

                 - m_dot [kg/s], the required flow rate

                 - if not explicitly set, PPI_current is AUTOMATICALLY RETRIEVED from:

                    "https://beta.bls.gov/dataQuery/find?st=0&r=20&s=popularity%3AD&fq=survey:
                    [pc]&more=0&q=Oil+and%20Gas%20Extraction"

                 - S_well, the expected success rate of the well drilling process (95% is the default value)

             The cost correlation is derived from:

                Adams et al, “Estimating the Geothermal Electricity Generation Potential of Sedimentary Basins Using
                genGEO (the generalizable GEOthermal techno-economic simulator)", ChemRxiv Prepr., 2021.

            Contingency cost increase (supposed 15%) -> 1.15
            Indirect cost increase (supposed 5%) -> 1.05

        """

        PPI_2002 = 84.0
        c_well, dc_CO2 = self.get_c_well()
        self.CAPEX_well = 1.15 * 1.05 * self.main_BHE.current_PPI / PPI_2002 * (c_well / self.main_BHE.s_well + dc_CO2)

    @abstractmethod
    def get_c_well(self):

        """

        This method must return the CAPEX of the well and the eventual sCO2 cost increase given the geometry and the
        design condition, as defined in the paper:

            Adams et al, “Estimating the Geothermal Electricity Generation Potential of Sedimentary Basins Using genGEO
            (the generalizable GEOthermal techno-economic simulator)", ChemRxiv Prepr., 2021.

        """

        return 0., 0.

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   OPTIMIZATION METHODS   ---------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    @abstractmethod
    def initial_guess(self):

        """
            This property must return the initial guesses for the parameter to be optimized
        """
        pass

    @property
    @abstractmethod
    def optimization_bounds(self):

        """
            This property must return the initial guesses for the parameter to be optimized
        """
        pass

    @abstractmethod
    def set_optimization_param(self, optimization_param):

        """
            This function must set the parameters provided by the optimizer
        """
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <---------------------   DISPLAY METHODS   ------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def plot_profiles(self):

        warnings.warn("\n\n!!IMPOSSIBLE TO DISPLAY THE PROFILE:\nThe heating section in use does not handle the desired function\n\n")
