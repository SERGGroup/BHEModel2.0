from main_code.constants import FREE_FEM_FOLDER
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import os, subprocess
import time

EDP_FILE_PATH = os.path.join(FREE_FEM_FOLDER, "Temperature Distribution.edp")
FREE_FEM_EXEC = "C:\Program Files (x86)\FreeFem++\FreeFem++.exe"



class FreeFEMOptions(ABC):

    def __init__(self):

        self.__init_calculation_folder(FREE_FEM_FOLDER)

    # FILE HANDLING METHODS

    def __init_calculation_folder(self, calculation_folder):

        self.__calculation_folder = self.check_and_create_dir(calculation_folder)

        self.__edp_file_path = os.path.join(self.__calculation_folder, "calculation_file.edp")
        self.__res_file_path = os.path.join(self.__calculation_folder, "result.txt")
        self.__plots_directory = self.check_and_create_dir(

            os.path.join(self.__calculation_folder, "plots")

        )

        self.other_calculation_folder_implementation()

    @staticmethod
    def check_and_create_dir(directory):

        if not os.path.isdir(directory):

            try:

                os.mkdir(directory)

            except:

                directory = FREE_FEM_FOLDER

        return directory

    @property
    def calculation_folder(self):

        return self.__calculation_folder

    @calculation_folder.setter
    def calculation_folder(self, new_calc_folder):

        self.__init_calculation_folder(new_calc_folder)

    @property
    def edp_file_path(self):

        return self.__edp_file_path

    @property
    def res_file_path(self):

        return self.__res_file_path

    def write_edp_file(self):

        with open(self.edp_file_path, "w") as edp_file:
            edp_file.write(self.edp_text())

    @abstractmethod
    def edp_text(self):

        pass

    @abstractmethod
    def other_calculation_folder_implementation(self):

        pass

class FreeFEMAnalyzer:

    def __init__(

            self, options: FreeFEMOptions,
            calculation_folder=FREE_FEM_FOLDER):

        self.options = options
        self.options.calculation_folder = calculation_folder

        self.reset_calculation_stats()

    def reset_calculation_stats(self):

        self.__ovr_calculation_time = 0.
        self.__ovr_calculation_calls = 0.

    def retrieve_result(self, clear_results=False):

        if os.path.isfile(self.options.res_file_path):

            with open(self.options.res_file_path, "r") as edp_file:

                lines = edp_file.readlines()

            if clear_results:
                os.remove(self.options.res_file_path)

            return lines

        else:

            return None

    def calculate(self, clear_results=False):

        self.options.write_edp_file()

        start_time = time.time()
        subprocess.call(

            '"{}" "{}"'.format(FREE_FEM_EXEC, self.options.edp_file_path),
            cwd=self.options.calculation_folder, shell=True,
            stdout=open(os.devnull, 'w'),
            stderr=subprocess.STDOUT

        )

        elapsed_time = time.time() - start_time
        self.__ovr_calculation_time += elapsed_time
        self.__ovr_calculation_calls += 1

        return self.retrieve_result(clear_results=clear_results)

    @property
    def mean_calculation_time(self):

        return self.__ovr_calculation_time / self.__ovr_calculation_calls