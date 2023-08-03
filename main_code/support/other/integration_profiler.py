

class IntegrationProfiler:

    def __init__(self, integrator):

        self.integrator_profile = list()
        self.integrator = integrator

    def step(self):

        if self.integrator.status == 'running':

            self.integrator.step()

            range_list = [self.integrator.t_old, self.integrator.t]
            range_list.sort()

            self.integrator_profile.append({

                "range": range_list,
                "dense_out": self.integrator.dense_output(),
                "error": (not self.integrator.status == 'failed')

            })

    def get_iteration_value(self, position):

        for step in self.integrator_profile:

            if step["range"][0] <= position <= step["range"][1]:

                return step["dense_out"](position)

    @property
    def t(self):
        return self.integrator.t

    @property
    def t_old(self):
        return self.integrator.t_old

    @property
    def status(self):
        return self.integrator.status

    @property
    def y(self):
        return self.integrator.y

    @property
    def y_old(self):
        return self.integrator.y_old
