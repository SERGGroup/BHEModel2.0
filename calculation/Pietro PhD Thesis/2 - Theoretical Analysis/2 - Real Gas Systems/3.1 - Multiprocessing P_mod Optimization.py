from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.power_plants.base_surface_power_plant import BaseSurfacePlant
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.constants import CALCULATION_FOLDER, os
from scipy.optimize import minimize
import concurrent.futures
import scipy.constants
import numpy as np
import math
import time
import gc

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "2 - Real Gas Systems",
    "0.support", "res"

)

support_folder = os.path.join(result_folder, "support")

base_filename = "res_{i}_{j}.npy"
reporting_points = [0.05, 0.25, 0.5, 0.75, 1]

t_rels = np.array([0.2, 1.1, 2])
n_temp = len(t_rels)

n_points = 5
grad_nd_nd_list = np.logspace(-1, 2, n_points)
dz_nd_list = np.logspace(-2, 2, n_points)
grad_nd_nd, dz_nd = np.meshgrid(grad_nd_nd_list, dz_nd_list)

fluid = "Methane"
bhe_in = PlantThermoPoint([fluid], [1])
bhe_in_base_si = PlantThermoPoint([fluid], [1], unit_system="MASS BASE SI")

bhe_in.set_variable("T", bhe_in.RPHandler.TC)
bhe_in.set_variable("P", bhe_in.RPHandler.PC)
rho_crit = bhe_in.get_variable("rho")
t_ambs = t_rels * bhe_in.RPHandler.TC


def update_values(in_state, in_base_state, t_amb, p_in, i, j, k):

    in_state.set_variable("T", t_amb)
    in_state.set_variable("P", p_in)
    in_state.copy_state_to(in_base_state)

    dpdt_in = in_base_state.get_derivative("P", "T", "rho")
    gamma_in = in_base_state.get_variable("CP/CV")
    v_in = 1 / in_base_state.get_variable("rho")
    cp_in = in_base_state.get_variable("cp")
    grad_nd_lim = cp_in * (1 - 1 / gamma_in) / (v_in * dpdt_in)
    grad_nd = grad_nd_nd[j, k] + grad_nd_lim

    grad = grad_nd * scipy.constants.g / cp_in
    depth = dz_nd[j, k] * (cp_in * in_state.get_variable("T")) / scipy.constants.g
    t_rocks = t_ambs[i] + grad * depth

    curr_bhe = SimplifiedBHE(in_state, dz_well=depth, t_rocks=t_rocks)
    new_base_surface = BaseSurfacePlant(curr_bhe)
    new_base_surface.prepare_calculation()
    new_base_surface.calculate()

    return new_base_surface


def evaluate_part(i, j):

    rep_pos = 0
    curr_value = i * n_points + j + 1
    max_digits = math.ceil(np.log10(n_points*n_temp))
    current_name = "{{:{}d}} of {{}}".format(max_digits).format(curr_value, n_points*n_temp)

    print("{} - Started".format(current_name))

    sub_start = time.time()

    sub_results = np.empty([2 * 3, n_points])
    sub_results[:] = np.nan

    tmp_bhe = bhe_in.duplicate()
    tmp_base_bhe = bhe_in_base_si.duplicate()

    tmp_bhe.set_variable("T", t_ambs[i])
    tmp_bhe.set_variable("rho", rho_crit)
    p_ref = tmp_bhe.get_variable("P")

    for k in range(n_points):

        def opt_func(x, opt_value=(0,)):

            p_curr = p_ref * x[0]
            base_surface = update_values(

                tmp_bhe, tmp_base_bhe,
                t_amb=t_ambs[i], p_in=p_curr,
                i=i, j=j, k=k

            )

            if opt_value[0] == 0:

                return -base_surface.w_dot_nd

            elif opt_value[0] == 1:

                return -base_surface.ex_dot_nd

            else:

                return -base_surface.eta_exs

        x0 = [2, 2, 0.8]
        for n in range(3):

            res = minimize(opt_func, np.array(x0[n]), args=[n])
            opt_res = -res.fun
            opt_res_x = res.x[0]

            res_test = -opt_func([0.9999], opt_value=(j,))
            if res_test > opt_res:
                opt_res_x = 0.9999
                opt_res = res_test

            res_test = -opt_func([1.00001], opt_value=(j,))
            if res_test > opt_res:
                opt_res_x = 1.00001
                opt_res = res_test

            sub_results[n, k] = opt_res_x
            sub_results[n + 3, k] = opt_res

        gc.collect()

        if k / n_points > reporting_points[rep_pos]:

            time_elapsed = time.time() - sub_start
            time_remaining = time_elapsed * (1 - reporting_points[rep_pos]) / reporting_points[rep_pos]

            print("{} - {:.0f}% done! -> {:.0f}s ({:.0f}min) elapsed, {:.0f}s ({:.0f}min) to-go".format(

                current_name, reporting_points[rep_pos]*100,
                time_elapsed, time_elapsed/60,
                time_remaining, time_remaining/60

            ))

            rep_pos += 1

    filename = os.path.join(support_folder, base_filename.format(i=i, j=j))
    np.save(filename, sub_results)

    time_elapsed = time.time() - sub_start

    print("{} - Done! -> {:.0f}s ({:.0f}min) elapsed".format(current_name, time_elapsed, time_elapsed / 60))

    return time_elapsed


if __name__ == '__main__':

    print("Calculation Started!!")

    with concurrent.futures.ProcessPoolExecutor() as executor:

        future_list = list()

        for i in range(n_temp):

            for j in range(n_points):

                time.sleep(0.5)
                future_list.append(executor.submit(evaluate_part, i, j))

        elapsed_time_list = list()
        for f in concurrent.futures.as_completed(future_list):

            try:
                elapsed_time_list.append(f.result())

            except:
                elapsed_time_list.append(np.nan)

    print("Calculation Completed!!")
    print("Saving Results ...")
    time.sleep(2)

    res_shape = np.append([2 * 3, n_temp], grad_nd_nd.shape)
    results = np.empty(res_shape)
    results[:] = np.nan

    for i in range(n_temp):

        for j in range(n_points):

            filename = os.path.join(support_folder, base_filename.format(i=i, j=j))

            if os.path.exists(filename):

                results[:, i, j, :] = np.load(filename)
                os.remove(filename)

    filename = os.path.join(result_folder, "opt_results.npy")
    np.save(filename, results)

    filename = os.path.join(result_folder, "depths_mesh.npy")
    np.save(filename, dz_nd)

    filename = os.path.join(result_folder, "grads_mesh.npy")
    np.save(filename, grad_nd_nd)

    filename = os.path.join(result_folder, "t_rels.npy")
    np.save(filename, t_rels)

    filename = os.path.join(result_folder, "elapsed_times.npy")
    np.save(filename, np.array(elapsed_time_list))

    print("Calculation Completed!!!")
