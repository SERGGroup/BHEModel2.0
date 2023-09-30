# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
import scipy.constants
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   INIT ARRAYS                         -------------------------------------> #
n_grad = 5
n_depth = 5
n_t_rel = 3

grad_nd_list = np.logspace(0, 2, n_grad + 1)[1:]
dz_nd_list = np.logspace(-2, 0, n_depth)
t_rel_list = np.logspace(np.log10(0.5), np.log10(5), n_t_rel)

grad_nd, dz_nd = np.meshgrid(grad_nd_list, dz_nd_list, indexing="ij")

base_shape = np.array(grad_nd.shape)
res_shape = np.append([len(t_rel_list)], base_shape)

grad_rocks = np.empty(res_shape)
depth = np.empty(res_shape)
t_rocks = np.empty(res_shape)

grad_rocks[:] = np.nan
depth[:] = np.nan
t_rocks[:] = np.nan

w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE                           -------------------------------------> #

fluid = PlantThermoPoint(["Methane"], [1])
v_rel_curr = 10**-2
p_in = fluid.

pbar = tqdm(desc="Calculating Points", total=len(t_rel_list) * len(grad_nd_list) * len(dz_nd_list))

for i in range(len(t_rel_list)):

    t_rel = t_rel_list[i]

    in_state = fluid.get_state(t=t_rel * fluid.t_crit, v=v_in)
    t_in = in_state.t

    grad_rocks[i, :, :] = grad_nd * scipy.constants.g / in_state.cp
    depth[i, :, :] = dz_nd * in_state.cp * t_in / scipy.constants.g
    t_rocks[i, :, :] = t_in + depth[i, :, :] * grad_rocks[i, :, :]

    for j in range(len(grad_nd_list)):

        for k in range(len(dz_nd_list)):

            states = evaluate_system(fluid, in_state, depth[i, j, k], t_rocks[i, j, k])
            surface_states = evaluate_surface(fluid, states, evaluate_with_flash=False)

            w_dot = (states[3].h - states[0].h)
            ex_dot = w_dot - states[0].t * (states[3].s - states[0].s)
            cf = 1 - t_in / t_rocks[i, j, k]

            w_dot_nds[i, j, k] = w_dot / (states[0].cp * states[0].t)
            ex_dot_nds[i, j, k] = ex_dot / (states[0].cp * states[0].t)
            eta_exs[i, j, k] = ex_dot / (w_dot * cf)

            w_dex_mins[i, j, k] = (surface_states[0].h - states[0].h) / w_dot
            w_dex_maxs[i, j, k] = (states[3].h - surface_states[1].h) / w_dot

            pbar.update(1)

pbar.close()