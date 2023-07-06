# %%------------
import scipy.integrate

from main_code.simplified_well.heating_sections import DefaultHeatingSection, AbstractHeatingSection
from main_code.support import PlantThermoPoint, retrieve_PPI
from main_code import constants

from EEETools.Tools.API.ExcelAPI.modules_importer import calculate_excel
from EEETools.Tools.API.Tools.main_tools import get_result_data_frames
from scipy.integrate import RK45
from scipy.optimize import Bounds
from scipy.constants import g
import numpy as np
from sty import ef
import warnings
import os
# %%------------
class SimplifiedBHE:

    def __init__(

            self, input_thermo_point: PlantThermoPoint, dz_well, T_rocks,
            heating_section=None, k_rocks=0.2, geo_flux=0.1, q_up_down=0.0,
            PPI=None, use_rk=False, d_inj=None, d_prod=None,
            discretize_p_losses=False

    ):

        # GEOTHERMAL FIELD CONDITION

        self.dz_well = dz_well      # unit: m .............. default: NONE (REQUIRED)
        self.T_rocks = T_rocks      # unit: °C ............. default: NONE (REQUIRED)
        self.k_rocks = k_rocks      # unit: W / (m K) ...... default: 0.2 W / (m K)
        self.geo_flux = geo_flux    # unit: W / m^2 ........ default: 0.1 W / m^2

        self.d_inj = d_inj          # unit: m .............. default: NONE (PRESSURE LOSSES IGNORED)
        self.d_prod = d_prod        # unit: m .............. default: NONE (PRESSURE LOSSES IGNORED)

        self.geo_gradient = None
        self.depth_optimization = False
        self.use_rk = use_rk
        self.discretize_p_losses = discretize_p_losses

        self.q_dot_up = - q_up_down
        self.q_dot_down = q_up_down

        self.c_well = 0.
        self.C0 = [0., 0.]
        self.P_loss = [0., 0.]

        self.__init_points(input_thermo_point)
        self.__init_heating_section(heating_section)
        self.__reset_control_elements(first_initialization=True)
        self.__ambient_condition = None

        self.__current_PPI = PPI

    def __init_points(self, input_thermo_point):

        self.points = [input_thermo_point]

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
