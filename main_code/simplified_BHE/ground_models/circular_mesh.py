from constants import GROUND_MESH_FOLDER
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import fortranformat as ff
import scipy.sparse.linalg
import scipy.sparse
import numpy as np
import os


class AbstractMeshGroundModel(ABC):

    def __init__(self, mesh_dim):

        self.mesh_points = np.zeros((2, 0, 0))
        self.boundary = dict()

        self.matrix = None
        self.base_vector = None
        self.vector = None
        self.T = None

        self.d_theta = 0.
        self.mesh_dim = mesh_dim

        self.__init_mesh()

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------   PROPERTIES   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def mesh_dim(self):
        return self.__mesh_dim

    @mesh_dim.setter
    def mesh_dim(self, mesh_dim_in):

        self.__mesh_dim = np.array(mesh_dim_in)
        self.check_mesh_dim()
        self.d_theta = 2 * np.pi / (self.n_x - 1)

    @property
    @abstractmethod
    def n_unknown(self):

        return 0.

    @property
    def n_x(self):
        return int(self.__mesh_dim[0])

    @property
    def n_y(self):
        return int(self.__mesh_dim[1])

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   MESH INITIALIZATION METHODS   ------------------------>
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def __init_mesh(self):

        mesh_file = self.__search_mesh_file()

        if mesh_file is not None:

            self.__open_mesh(mesh_file)

        else:

            self.generate_mesh()
            self.__store_mesh(store_in_mesh_repository=True)

    @abstractmethod
    def generate_mesh(self):

        pass

    @abstractmethod
    def init_boundary(self):

        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <--------------   MATRIX INITIALIZATION & SOLUTION   --------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @abstractmethod
    def solve(self, DT, grad=1):

        pass

    def evaluate_standard_solution(self):

        if self.matrix is None:
            self.__init_matrix()

        self.init_vector()
        sol = scipy.sparse.linalg.spsolve(self.matrix.tocsr(copy=True), self.vector.tocsr(copy=True))

        self.append_solution(sol)
        self.additional_calculation()

    def __init_matrix(self):

        n_unknown = self.n_unknown
        self.matrix = scipy.sparse.lil_matrix((n_unknown, n_unknown))
        self.base_vector = scipy.sparse.lil_matrix((n_unknown, 1))

        # initialize diagonal
        for i in range(self.n_x):
            self.matrix[i, i] = 0.

        # append values
        for i in range(self.n_x - 1):

            for j in range(self.n_y - 1):

                self.init_matrix_row(i, j)

    @abstractmethod
    def init_matrix_row(self, i, j):

        pass

    @abstractmethod
    def init_vector(self):

        pass

    @abstractmethod
    def append_solution(self, sol):

        pass

    @abstractmethod
    def additional_calculation(self):

        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   MESH STORE & RETRIEVE   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def __search_mesh_file(self):

        for filename in os.listdir(os.path.join(GROUND_MESH_FOLDER, self.mesh_dir_name)):

            file = os.path.join(GROUND_MESH_FOLDER, self.mesh_dir_name, filename)

            if os.path.isfile(file):

                if self.is_usable_file(filename):
                    return file

        return None

    @abstractmethod
    def is_usable_file(self, filename):

        pass

    @property
    @abstractmethod
    def filename(self):

        return ""

    @property
    @abstractmethod
    def mesh_dir_name(self):

        return ""

    def __open_mesh(self, mesh_file):

        counter = 0
        n_x = self.n_x
        n_y = self.n_y

        self.mesh_points = np.zeros((2, n_x, n_y))
        self.init_boundary()

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

                self.mesh_points[:, i, j] = self.modify_points(np.float64(data), reverse_modification=True)
                counter += 1

    def __store_mesh(self, store_in_mesh_repository=False):

        filename = self.filename

        if store_in_mesh_repository:

            file_path = os.path.join(GROUND_MESH_FOLDER, self.mesh_dir_name, filename)

        else:

            CURRENT_FOLDER = os.path.dirname(os.path.abspath(__file__))
            file_path = os.path.join(CURRENT_FOLDER, filename)

        with open(file_path, 'w') as f:

            f.write(str(self))

    @staticmethod
    def __read_mesh_file(mesh_file):

        points = None

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

                if points is not None:

                    points = np.vstack([points, np.float64(data)])

                else:

                    points = np.float64(data)

        return points

    @abstractmethod
    def modify_points(self, value, reverse_modification=False):
        pass

    def __str__(self):

        self_str = ""
        line_format = ff.FortranRecordWriter('(E16.8,A,E16.8)')

        for i in range(int(self.__mesh_dim[0])):

            for j in range(int(self.__mesh_dim[1])):

                modified_point = self.modify_points(self.mesh_points[:, i, j])

                self_str += "{}\n".format(

                    line_format.write([

                        modified_point[0], ";",
                        modified_point[1]

                    ]),

                )

        return self_str

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------------   SUPPORT METHODS   ------------------------------>
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @abstractmethod
    def check_mesh_dim(self):

        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   PLOT & EXPORT METHODS   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def plot_mesh(self):

        for i in range(self.n_x):
            plt.plot(self.mesh_points[0, i, :], self.mesh_points[1, i, :], '#929591', linewidth=0.5)

        for j in range(self.n_y):
            plt.plot(self.mesh_points[0, :, j], self.mesh_points[1, :, j], '#808080', linewidth=0.5)

        for key in self.boundary.keys():

            points = self.boundary[key]
            plt.plot(points[0, :], points[1, :], '#DBB40C')

        plt.show()

    def plot_field(self):

        if not self.T is None:

            plt.contourf(self.mesh_points[0], self.mesh_points[1], self.T, 20, cmap='RdGy_r')
            plt.colorbar()
            plt.show()


class CircularMeshGroundModel(AbstractMeshGroundModel):

    def __init__(self, mesh_dim, d_tube, d_ratio, alpha, use_square_adaptor=False):

        self.r_tube = d_tube / 2
        self.d_ratio = d_ratio
        self.max_pos = self.r_tube * d_ratio
        self.alpha = alpha

        self.DT = 0.
        self.grad = 0.

        super().__init__(mesh_dim)
        self.__init_square_adaptor(use_square_adaptor)

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------   PROPERTIES   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def n_unknown(self):

        n_nodes = self.n_x * self.n_y
        n_boundary = 2 * self.n_x + self.n_y - 2
        return n_nodes - n_boundary

    @property
    def filename(self):

        return '{}-{}-{}.dat'.format(

            self.mesh_dim[0],
            self.mesh_dim[1],
            self.alpha

        )

    @property
    def mesh_dir_name(self):

        return "circ_mesh_repo"

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   MESH INITIALIZATION METHODS   ------------------------>
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def generate_mesh(self):

        self.init_boundary()
        self.__init_other_points()

    def init_boundary(self):

        n_x = self.n_x
        n_y = self.n_y

        self.boundary = {

            "Top": np.zeros((2, n_x)),
            "Bottom": np.zeros((2, n_x)),
            "Left": np.zeros((2, n_y)),
            "Right": np.zeros((2, n_y))

        }

        for i in range(n_x):
            self.boundary["Top"][:, i] = self.__init_top_boundary(i / (n_x - 1))
            self.boundary["Bottom"][:, i] = self.__init_bottom_boundary(i / (n_x - 1))

        for j in range(n_y):
            point = self.__init_side_boundary(j / (n_y - 1))
            self.boundary["Left"][:, j] = point
            self.boundary["Right"][:, j] = point

    def __init_other_points(self):

        n_x = self.n_x
        n_y = self.n_y
        self.mesh_points = np.zeros((2, n_x, n_y))

        for i in range(n_x):

            top_point = self.boundary["Top"][:, i]
            bottom_point = self.boundary["Bottom"][:, i]

            perc = 0
            dperc = (1 - self.alpha) / (1 - pow(self.alpha, n_y - 1))

            for j in range(n_y):
                self.mesh_points[:, i, j] = bottom_point * (1 - perc) + top_point * perc

                perc = perc + dperc
                dperc = self.alpha * dperc

    def __init_top_boundary(self, perc):

        theta = 2 * np.pi * perc
        point = np.array((np.cos(theta), np.sin(theta)))

        return point * self.max_pos

    def __init_bottom_boundary(self, perc):

        theta = 2 * np.pi * perc
        point = np.array((np.cos(theta), np.sin(theta)))

        return point * self.r_tube

    def __init_side_boundary(self, perc):

        return np.array((self.r_tube * (1 - perc) + self.max_pos * perc, 0))

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <--------------   MATRIX INITIALIZATION & SOLUTION   --------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def solve(self, DT, grad=1):

        self.DT = DT
        self.grad = grad

        if self.square_adaptor is None:

            self.evaluate_standard_solution()

        else:

            self.__evaluate_square_adaptor_solution()

    def __get_index(self, i, j):

        i = np.mod(i, self.n_x - 1)

        if j == self.n_y - 1:

            return -1

        else:

            j = j - 1
            return i * (self.n_y - 2) + j

    def init_matrix_row(self, i, j):

        i = np.mod(i, self.n_x - 1)
        n_node = self.__get_index(i, j)
        r_curr = self.mesh_points[0][0][j]
        r_up = self.mesh_points[0][0][j + 1]

        # upper node (j = j + 1)
        n_up = self.__get_index(i, j + 1)
        g_up = self.d_theta / np.log(r_up / r_curr)

        if j == 0:

            self.base_vector[n_up] = -g_up
            self.matrix[n_up, n_up] -= g_up

        elif j == self.n_y - 2:

            self.base_vector[n_node] = -g_up
            self.matrix[n_node, n_node] -= g_up

        else:

            # left node (i = i + 1)
            n_left = self.__get_index(i + 1, j)
            r_down = self.mesh_points[0][0][j - 1]
            g_left = (r_up - r_down) / (2 * r_curr * self.d_theta)

            # As the matrix must be symmetric the calculation are performed only once and the calculated resistance
            # and written twice in the matrix

            # append off-diagonal values
            self.matrix[n_node, n_up] = g_up
            self.matrix[n_up, n_node] = g_up
            self.matrix[n_node, n_left] = g_left
            self.matrix[n_left, n_node] = g_left

            # update diagonal values
            self.matrix[n_up, n_up] -= g_up
            self.matrix[n_node, n_node] -= g_up
            self.matrix[n_node, n_node] -= g_left
            self.matrix[n_left, n_left] -= g_left

    def init_vector(self):

        self.T = np.zeros((self.n_x, self.n_y))
        self.vector = scipy.sparse.lil_matrix((self.n_unknown, 1))

        for i in range(self.n_x):
            node_up = self.__get_index(i, self.n_y - 2)
            node_down = self.__get_index(i, 1)

            perc = i / (self.n_x - 1)
            dT_down = - self.DT
            dT_up = - np.sin(2 * np.pi * perc) * self.grad * self.max_pos

            self.T[i, 0] = dT_down
            self.T[i, -1] = dT_up

            self.vector[node_down] = self.base_vector[node_down] * dT_down
            self.vector[node_up] = self.base_vector[node_up] * dT_up

    def append_solution(self, sol):

        for i in range(self.n_x):

            for j in range(1, self.n_y - 1):

                n_pos = self.__get_index(i, j)
                if not n_pos == -1:
                    self.T[i, j] = sol[n_pos]

    def additional_calculation(self):

        self.__evaluate_internal_gradient(calc_quadratic=True)
        self.__evaluate_internal_flux()

    def __evaluate_internal_gradient(self, calc_quadratic=False):

        dr_int = self.mesh_points[0][0][1] - self.mesh_points[0][0][0]

        if calc_quadratic:

            dr_ext = self.mesh_points[0][0][2] - self.mesh_points[0][0][0]

            __T_1 = self.T[0:-1, 0]
            __T_2 = self.T[0:-1, 1]
            __T_3 = self.T[0:-1, 2]

            __dT_3_dot = (__T_3 - __T_1) / dr_ext
            __dT_2_dot = (__T_2 - __T_1) / dr_int

            __a = (__dT_3_dot - __dT_2_dot) / (dr_ext - dr_int)

            gradients = __dT_3_dot - __a * dr_ext

        else:

            gradients = (self.T[0:-1, 1] - self.T[0:-1, 0]) / dr_int

        self.int_grad = abs(np.mean(gradients))
        self.gamma = self.int_grad / self.grad

    def __evaluate_internal_flux(self):

        self.flux = 0.
        counter = 0

        # append values
        for i in range(self.n_x - 1):

            n = self.__get_index(i, 1)
            dT = self.T[i, 1] - self.T[i, 0]
            self.flux -= self.base_vector[n, 0] * dT
            counter += 1

        self.q_ratio_ND = self.flux / (2 * np.pi * self.r_tube * self.grad)
        self.q_ratio = self.flux / (counter * self.grad)

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------   SQUARE ADAPTOR SOLUTION   -------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def __init_square_adaptor(self, use_square_adaptor):

        self.square_adaptor = None

        if use_square_adaptor:
            self.square_adaptor = SquareAdaptor(self)

            # <------------------------------------------------------------------------->

    def __evaluate_square_adaptor_solution(self):

        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   MESH STORE & RETRIEVE   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def modify_points(self, points, reverse_modification=False):

        if reverse_modification:

            return ((self.d_ratio - 1) * points + 1) * self.r_tube

        else:

            return (points / self.r_tube - 1) / (self.d_ratio - 1)

    def is_usable_file(self, filename):

        try:

            data_elements = filename.strip(".dat")

            if "-" in data_elements:

                data_elements = data_elements.split("-")
                dim_check = [int(data_elements[0]), int(data_elements[1])] == self.__mesh_dim

                if len(data_elements) == 4:

                    ratio_check = float(data_elements[2]) == float(self.d_ratio)
                    alpha = float(data_elements[2])

                else:

                    ratio_check = True
                    alpha = float(data_elements[3])

                alpha_check = (alpha == self.alpha)

                return dim_check[0] and dim_check[1] and ratio_check and alpha_check

        except:

            return False

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------------   SUPPORT METHODS   ------------------------------>
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def check_mesh_dim(self):

        if not (self.mesh_dim[0] - 1) % 8 == 0:
            self.mesh_dim[0] = np.ceil(self.mesh_dim[0] / 8) * 8 + 1

        for i in range(2):
            if self.mesh_dim[i] < 4:
                self.mesh_dim[i] = 4


class SquareAdaptor(AbstractMeshGroundModel):

    def __init__(self, circular_mesh):

        self.circular_mesh = circular_mesh
        super().__init__(self.__evaluate_mesh_dim())

    def __evaluate_mesh_dim(self):
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------   PROPERTIES   --------------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    @property
    def n_unknown(self):
        return 0.

    @property
    def filename(self):
        return ""

    @property
    def mesh_dir_name(self):
        return os.path.join("circ_mesh_repo", "square_adaptor")

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------   MESH INITIALIZATION METHODS   ------------------------>
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def generate_mesh(self):
        pass

    def init_boundary(self):
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <--------------   MATRIX INITIALIZATION & SOLUTION   --------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def solve(self, DT, grad=1):
        pass

    def init_matrix_row(self, i, j):
        pass

    def init_vector(self):
        pass

    def append_solution(self, sol):
        pass

    def additional_calculation(self):
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <-------------------   MESH STORE & RETRIEVE   --------------------------->
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def is_usable_file(self, filename):
        pass

    def modify_points(self, value, reverse_modification=False):
        pass

    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->
    # <----------------------   SUPPORT METHODS   ------------------------------>
    # <------------------------------------------------------------------------->
    # <------------------------------------------------------------------------->

    def check_mesh_dim(self):
        pass


def mesh_test():

    mesh = CircularMeshGroundModel(d_tube=2, d_ratio=500, mesh_dim=[200, 400], alpha=1.01)
    mesh.solve(100)
    print(mesh.flux)


def d_ratio_analysis():

    x = list()
    y = list()
    y_1 = list()
    y_2 = list()

    div = 1

    for i in range(3):

        d_ratio = np.power(10, (i + 1) * div + 2)

        mesh = CircularMeshGroundModel(d_tube=2, d_ratio=d_ratio, mesh_dim=[200, 400], alpha=1.01)
        mesh_1 = CircularMeshGroundModel(d_tube=2, d_ratio=d_ratio, mesh_dim=[100, 200], alpha=1.04)
        mesh_2 = CircularMeshGroundModel(d_tube=2, d_ratio=d_ratio, mesh_dim=[50, 100], alpha=1.06)

        mesh.solve(100)
        mesh_1.solve(100)
        mesh_2.solve(100)

        print("{}->{} - {} - {}".format(d_ratio, mesh.flux, mesh_1.flux, mesh_2.flux))

        x.append(d_ratio)
        y.append(mesh.flux)
        y_1.append(mesh_1.flux)
        y_2.append(mesh_2.flux)

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.plot(x, y, label="mesh_0")
    ax.plot(x, y_1, label="mesh_1")
    ax.plot(x, y_2, label="mesh_2")

    ax.set_xscale('log')
    plt.show()


def dt_analysis():

    fig = plt.figure()
    ax = fig.add_subplot()

    for d_tube in range(2, 10, 2):

        x = list()
        y = list()
        mesh = CircularMeshGroundModel(d_tube=d_tube, d_ratio=1000, mesh_dim=[50, 100], alpha=1.06)

        for dT in range(5, 150):

            mesh.solve(dT)
            print("{} -> {}".format(dT, mesh.q_ratio_ND))

            x.append(dT)
            y.append(mesh.q_ratio_ND)

        ax.plot(x, y, label=str(d_tube))

    plt.show()


if __name__ == "__main__":

    dt_analysis()