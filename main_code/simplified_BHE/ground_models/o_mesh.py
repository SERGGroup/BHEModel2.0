from constants import GROUND_MESH_FOLDER
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
import fortranformat as ff
import numpy as np
import subprocess
import os


class GroundModel:
    pass


class GroundMesh:

    def __init__(self, d_tube, d_ratio, mesh_dim, max_iteration=8000, iteration_toll=10 ** -4,
                 display_mesh_iteration=False):

        self.r_tube = d_tube / 2
        self.d_ratio = d_ratio
        self.max_pos = self.r_tube * d_ratio

        self.mesh_generation_options = {

            "max_iteration": max_iteration,
            "iteration_toll": iteration_toll,
            "display_mesh_iteration": display_mesh_iteration

        }

        self.__append_mesh_dim(mesh_dim)
        self.__init_mesh()

    def __init_mesh(self):

        mesh_file = self.__search_mesh_file()

        if mesh_file is not None:

            self.__open_mesh(mesh_file)

        else:

            self.__generate_mesh()

    def __search_mesh_file(self):

        for filename in os.listdir(os.path.join(GROUND_MESH_FOLDER, "o_mesh_repo")):

            file = os.path.join(GROUND_MESH_FOLDER, "o_mesh_repo", filename)

            if os.path.isfile(file):

                if self.__is_usable_file(filename):
                    return file

        return None

    def __is_usable_file(self, filename):

        try:

            data_elements = filename.strip(".dat")

            if "-" in data_elements:
                data_elements = data_elements.split("-")
                dim_check = [int(data_elements[0]), int(data_elements[1])] == self.mesh_dim
                ratio_check = float(data_elements[2]) == float(self.d_ratio)

                return dim_check[0] and dim_check[1] and ratio_check

        except:

            return False

    def __open_mesh(self, mesh_file):

        counter = 0
        n_x = int(self.mesh_dim[0])
        n_y = int(self.mesh_dim[1])
        self.mesh_points = np.zeros((2, n_x, n_y))
        self.__init_boundary()

        with open(mesh_file, 'r') as file:

            for line in file.readlines():

                data = line.strip("\n")

                if "\t" in data:

                    data = line.strip()
                    data = data.split("\t")

                elif ";" in data:

                    data = line.strip()
                    data = data.split(";")

                else:

                    data = data.split(" ")

                i = int(np.floor(counter / n_y))
                j = int(np.mod(counter, n_y))
                self.mesh_points[:, i, j] = np.float64(data) * self.r_tube
                counter += 1

        self.plot()

    def __store_mesh(self, store_in_mesh_repository=False):

        filename = '{}-{}-{:.3E}.dat'.format(

            self.mesh_dim[0],
            self.mesh_dim[1],
            self.d_ratio

        )

        if store_in_mesh_repository:

            file_path = os.path.join(GROUND_MESH_FOLDER, "o_mesh_repo", filename)

        else:

            CURRENT_FOLDER = os.path.dirname(os.path.abspath(__file__))
            file_path = os.path.join(CURRENT_FOLDER,filename)

        with open(file_path, 'w') as f:

            f.write(str(self))

        self.plot()

        return filename

    def __generate_mesh(self, generate_in_python=False):

        self.__init_boundary()
        self.__init_empty_mesh()

        if generate_in_python:

            self.__smooth_mesh(

                self.mesh_generation_options["max_iteration"],
                self.mesh_generation_options["iteration_toll"],
                self.mesh_generation_options["display_mesh_iteration"]

            )

            self.__store_mesh(store_in_mesh_repository=True)

        else:

            # Fortran Smoothing
            filename = self.__store_mesh()
            exe_path = os.path.join(GROUND_MESH_FOLDER, "mesh_smoothing.exe")
            CURRENT_FOLDER = os.path.dirname(os.path.abspath(__file__))

            args = "{} {}".format(exe_path, filename)
            subprocess.call(args)
            os.replace(os.path.join(CURRENT_FOLDER, filename), os.path.join(GROUND_MESH_FOLDER, "o_mesh_repo", filename))

            self.__open_mesh(os.path.join(GROUND_MESH_FOLDER, "o_mesh_repo", filename))

    def __append_mesh_dim(self, mesh_dim):

        self.mesh_dim = np.array(mesh_dim)
        self.check_mesh_dim()

    def __init_boundary(self):

        n_x = int(self.mesh_dim[0])
        n_y = int(self.mesh_dim[1])

        self.boundary = {

            "Top": np.zeros((2, n_x)),
            "Bottom": np.zeros((2, n_x)),
            "Left": np.zeros((2, n_y)),
            "Right": np.zeros((2, n_y))

        }

        for i in range(n_x):
            self.boundary["Top"][:, i] = self.calc_top_boundary(i / (n_x - 1))
            self.boundary["Bottom"][:, i] = self.calc_bottom_boundary(i / (n_x - 1))

        for j in range(n_y):
            point = self.calc_side_boundary(j / (n_y - 1))
            self.boundary["Left"][:, j] = point
            self.boundary["Right"][:, j] = point

    def __init_empty_mesh(self):

        n_x = int(self.mesh_dim[0])
        n_y = int(self.mesh_dim[1])
        self.mesh_points = np.zeros((2, n_x, n_y))

        for i in range(n_x):

            top_point = self.boundary["Top"][:, i]
            bottom_point = self.boundary["Bottom"][:, i]

            for j in range(n_y):

                perc = j / (n_y - 1)
                self.mesh_points[:, i, j] = bottom_point * (1 - perc) + top_point * perc

    def __smooth_mesh(self, max_iteration, toll, display):

        pbar = tqdm(total=100, ncols=200, unit="%")

        self.conv_history = np.zeros((2, max_iteration))
        min_conv = 1

        for t in range(max_iteration):

            self.__update_point_smoothness(t)

            min_conv = np.min([min_conv, np.max(self.conv_history[:, t])])

            if display:
                self.plot()

            if min_conv < toll:
                break

            toll_perc = int(np.log(min_conv) / np.log(toll) * 100)
            d_perc = toll_perc - pbar.n
            pbar.set_description("Smoothing Mesh (Iter Limit: {:5.2f}%)".format(t / max_iteration * 100))
            pbar.update(d_perc)

        pbar.close()

    def __update_point_smoothness(self, t):

        n_x = int(self.mesh_dim[0])
        n_y = int(self.mesh_dim[1])
        new_points = np.zeros((2, n_x, n_y))

        x = self.mesh_points[0, :, :]
        y = self.mesh_points[1, :, :]

        # Boundary Conditions
        new_points[:, :, 0] = self.mesh_points[:, :, 0]
        new_points[:, :, n_y - 1] = self.mesh_points[:, :, n_y - 1]

        new_points[1, 0, :] = self.mesh_points[1, 0, :]
        new_points[1, n_x - 1, :] = self.mesh_points[1, n_x - 1, :]

        for i in range(1, n_x - 1):

            for j in range(1, n_y - 1):

                a = 1 / 4 * ((x[i, j + 1] - x[i, j - 1]) ** 2 + (y[i, j + 1] - y[i, j - 1]) ** 2)
                b = 1 / 16 * ((x[i + 1, j] - x[i - 1, j]) * (x[i, j + 1] - x[i, j - 1]) + (y[i + 1, j] - y[i - 1, j]) * (y[i, j + 1] - y[i, j - 1]))
                g = 1 / 4 * ((x[i + 1, j] - x[i - 1, j]) ** 2 + (y[i + 1, j] - y[i - 1, j]) ** 2)

                smoothing = 1 / (2 * (a + g + np.float_power(10, -9)))

                for n in range(2):

                    dn_tot = self.mesh_points[n, i + 1, j + 1] - self.mesh_points[n, i - 1, j + 1] - self.mesh_points[n, i + 1, j - 1] + self.mesh_points[n, i - 1, j - 1]
                    dn_x = self.mesh_points[n, i + 1, j] + self.mesh_points[n, i - 1, j]
                    dn_y = self.mesh_points[n, i, j + 1] + self.mesh_points[n, i, j - 1]

                    new_points[n, i, j] = - smoothing * (2 * b * dn_tot - a * dn_x - g * dn_y)
                    variation = np.abs(new_points[n, i, j] - self.mesh_points[n, i, j])
                    self.conv_history[n, t] = max(self.conv_history[n, t], variation)

        x_mean = (new_points[0, n_x - 2, :] + new_points[0, 1, :]) / 2
        new_points[0, n_x - 1, :] = x_mean
        new_points[0, n_x - 2, :] = x_mean
        new_points[0, 0, :] = x_mean
        new_points[0, 1, :] = x_mean

        self.mesh_points = new_points

    def plot(self):

        for i in range(int(self.mesh_dim[0])):
            plt.plot(self.mesh_points[0, i, :], self.mesh_points[1, i, :], '#929591', linewidth=0.5)

        for j in range(int(self.mesh_dim[1])):
            plt.plot(self.mesh_points[0, :, j], self.mesh_points[1, :, j], '#808080', linewidth=0.5)

        for key in self.boundary.keys():
            points = self.boundary[key]
            plt.plot(points[0, :], points[1, :], '#DBB40C')

        plt.show()

    def check_mesh_dim(self):

        if not (self.mesh_dim[0] - 1) % 8 == 0:

            self.mesh_dim[0] = np.ceil(self.mesh_dim[0] / 8) * 8 + 1

        for i in range(2):

            if self.mesh_dim[i] < 4:

                self.mesh_dim[i] = 4

    def calc_top_boundary(self, perc):

        point = np.zeros(2)

        if perc < 1 / 8:

            alpha = perc * 8
            point[:] = (self.max_pos, self.max_pos * alpha)

        elif perc < 3 / 8:

            alpha = (perc - 1 / 8) * 4
            point[:] = (self.max_pos * (1 - 2 * alpha), self.max_pos)

        elif perc < 5 / 8:

            alpha = (perc - 3 / 8) * 4
            point[:] = (- self.max_pos, self.max_pos * (1 - 2 * alpha))

        elif perc < 7 / 8:

            alpha = (perc - 5 / 8) * 4
            point[:] = (self.max_pos * (2 * alpha - 1), - self.max_pos)

        else:

            alpha = (perc - 7 / 8) * 8
            point[:] = (self.max_pos, self.max_pos * (alpha - 1))

        return point

    def calc_bottom_boundary(self, perc):

        theta = 2 * np.pi * perc
        point = np.array((np.cos(theta), np.sin(theta)))

        return point * self.r_tube

    def calc_side_boundary(self, perc):

        return np.array((self.r_tube * (1 - perc) + self.max_pos * perc, 0))

    def __str__(self):

        self_str = ""
        line_format = ff.FortranRecordWriter('(E16.8,A,E16.8)')

        for i in range(int(self.mesh_dim[0])):

            for j in range(int(self.mesh_dim[1])):

                self_str += "{}\n".format(

                    line_format.write([

                        self.mesh_points[0, i, j] / self.r_tube, ";",
                        self.mesh_points[1, i, j] / self.r_tube

                    ]),

                )

        return self_str


if __name__ == "__main__":

    mesh = GroundMesh(d_tube=2, d_ratio=100, mesh_dim=[400, 200])

    #for filename in os.listdir(GROUND_MESH_FOLDER):

     #   file = os.path.join(GROUND_MESH_FOLDER, filename)

      #  if os.path.isfile(file):

        #    data_elements = filename.strip(".dat")

         #   if "-" in data_elements:

          #      data_elements = data_elements.split("-")
           #     ground_mesh = GroundMesh(d_tube=2, d_ratio=float(data_elements[2]), __mesh_dim=[int(data_elements[0]), int(data_elements[1])])
