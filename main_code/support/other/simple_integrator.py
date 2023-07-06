from scipy.interpolate import interp1d
from copy import deepcopy
import numpy as np

class SimpleIntegrator:

    def __init__(self, integration_funct, pos_0, x0, pos_end, n_steps):

        self.integration_funct = integration_funct
        self.pos_0 = pos_0
        self.pos_end = pos_end
        self.x0 = x0
        self.n_steps = n_steps

        self.dpos =  (pos_end - pos_0) / (n_steps - 1)
        self.__curr_step = 0
        self.__curr_x = np.array(x0)
        self.__old_x = np.array(x0)

    @property
    def t(self):
        return self.pos_0 + self.dpos * self.__curr_step

    @property
    def t_old(self):
        return self.pos_0 + self.dpos * (self.__curr_step - 1)

    @property
    def status(self):

        if self.__curr_step < self.n_steps:
            return 'running'

        return 'finished'

    @property
    def y(self):

        return self.__curr_x

    @property
    def y_old(self):

        return self.__old_x

    def step(self):

        if self.__curr_step < self.n_steps:

            dx = self.integration_funct(self.t, self.y)

            self.__old_x = deepcopy(self.__curr_x)
            self.__curr_x += np.array(dx) * self.dpos
            self.__curr_step += 1

    def dense_output(self):

        x = np.array([self.t_old, self.t])
        y = np.stack((np.stack(self.y_old), np.array(self.y))).T
        return interp1d(x, y)