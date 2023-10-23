from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.constants import CALCULATION_FOLDER, os
import concurrent.futures
import scipy.constants
import numpy as np
import math
import time
import gc

n_grad = 150
n_depth = 3
n_t_rel = 40
n_p_rel = 40

p_rel_list = np.logspace(0.5, 2, n_p_rel)
t_rel_list = np.linspace(0.5, 2, n_t_rel)

# p_rel_list = [1, 10**1, 10**2]
# t_rel_list = [0.5, 0.75, 1.1, 2]

dz_nd_list = np.logspace(-2, 0, n_depth)
grad_nd_nd_list = np.logspace(-2, 4.5, n_grad)
grad_nd_nd, dz_nd = np.meshgrid(grad_nd_nd_list, dz_nd_list, indexing="ij")

n_processes = len(p_rel_list) * len(t_rel_list)
bhe_in = PlantThermoPoint(["Methane"], [1])
in_state = PlantThermoPoint(["Methane"], [1], unit_system="MASS BASE SI")

states = list()
for i in range(4):
    states.append(in_state.duplicate())

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "2 - Real Gas Systems",
    "0.support", "res"

)

support_folder = os.path.join(result_folder, "support")

base_filename = "res_{key}_{a}_{i}.npy"
reporting_points = [0.05, 0.25, 0.5, 0.75, 1]


def init_result_dict(res_shape):

    result_dict = {

        "grad_nd_rocks": np.empty(res_shape),
        "grad_rocks": np.empty(res_shape),
        "depth": np.empty(res_shape),
        "alpha_in": np.empty(res_shape),
        "t_rocks": np.empty(res_shape),
        "w_dot_nds": np.empty(res_shape),
        "ex_dot_nds": np.empty(res_shape),
        "eta_exs": np.empty(res_shape),
        "w_dex_mins": np.empty(res_shape),
        "w_dex_maxs": np.empty(res_shape)

    }

    for key in result_dict.keys():
        result_dict[key][:] = np.nan

    return result_dict


def evaluate_points(a_curr, i_curr):

    rep_pos = 0
    curr_value = a_curr * len(t_rel_list) + i_curr + 1
    max_digits = math.ceil(np.log10(n_processes))
    current_name = "{{:{}d}} of {{}}".format(max_digits).format(curr_value, n_processes)

    print("{} - Started".format(current_name))
    sub_start = time.time()

    res_shape = np.array([len(grad_nd_nd_list), len(dz_nd_list)])
    result_dict = init_result_dict(res_shape)

    p_rel = p_rel_list[a_curr]
    t_in = t_rel_list[i_curr] * in_state.RPHandler.TC
    p_in = p_rel * in_state.RPHandler.PC
    in_state.set_variable("T", t_in)
    in_state.set_variable("P", p_in)
    in_state.copy_state_to(bhe_in)

    cp_in = in_state.get_variable("cp")
    result_dict["depth"][:, :] = dz_nd * in_state.get_variable("cp") * t_in / scipy.constants.g

    for k in range(len(dz_nd_list)):

        bhe = SimplifiedBHE(bhe_in, dz_well=result_dict["depth"][0, k], t_rocks=50)
        bhe.update()

        t_input = bhe.points[0].get_variable("T")
        t_down = bhe.points[1].get_variable("T")

        alpha_in = cp_in * (t_down - t_input) / (scipy.constants.g * result_dict["depth"][0, k])
        result_dict["alpha_in"][:, k] = alpha_in
        result_dict["grad_nd_rocks"][:, k] = grad_nd_nd[:, k] + alpha_in

        grad_rock = result_dict["grad_nd_rocks"][:, k] * scipy.constants.g / in_state.get_variable("cp")
        t_rock = t_in + result_dict["depth"][:, k] * grad_rock
        result_dict["grad_rocks"][:, k] = grad_rock
        result_dict["t_rocks"][:, k] = t_rock

        bhe_subsub_list = list()
        for j in range(len(grad_nd_nd_list)):

            bhe = SimplifiedBHE(

                bhe_in, dz_well=result_dict["depth"][j, k],
                t_rocks=result_dict["t_rocks"][j, k] - 273.15

            )
            bhe.update()
            bhe_subsub_list.append(bhe)

            for n in range(4):
                bhe.points[n].copy_state_to(states[n])

            surface_states = list()
            surface_states.append(states[3].duplicate())
            surface_states.append(states[3].duplicate())

            surface_states[0].set_variable("P", states[3].get_variable("P"))
            surface_states[0].set_variable("S", states[0].get_variable("S"))

            surface_states[1].set_variable("P", states[0].get_variable("P"))
            surface_states[1].set_variable("S", states[3].get_variable("S"))

            w_dot = states[3].get_variable("H") - states[0].get_variable("H")
            ex_dot = w_dot - states[0].get_variable("T") * (states[3].get_variable("S") - states[0].get_variable("S"))

            if (not (w_dot < 0 or ex_dot < 0)) and (not states[-1].get_variable("H") < -1e6):

                cf = 1 - states[0].get_variable("T") / states[2].get_variable("T")
                w_dot_nd = w_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))
                ex_dot_nd = ex_dot / (states[0].get_variable("CP") * states[0].get_variable("T"))

                result_dict["w_dot_nds"][j, k] = w_dot_nd
                result_dict["ex_dot_nds"][j, k] = ex_dot_nd
                result_dict["eta_exs"][j, k] = ex_dot / (w_dot * cf)

                w_dex_min = (surface_states[1].get_variable("H") - states[0].get_variable("H")) / w_dot
                w_dex_max = (states[3].get_variable("H") - surface_states[1].get_variable("H")) / w_dot
                result_dict["w_dex_mins"][j, k] = w_dex_min
                result_dict["w_dex_maxs"][j, k] = w_dex_max

            gc.collect()

        if k / len(dz_nd_list) > reporting_points[rep_pos]:

            time_elapsed = time.time() - sub_start
            time_remaining = time_elapsed * (1 - reporting_points[rep_pos]) / reporting_points[rep_pos]

            print("{} - {:.0f}% done! -> {:.0f}s ({:.0f}min) elapsed, {:.0f}s ({:.0f}min) to-go".format(

                current_name, reporting_points[rep_pos]*100,
                time_elapsed, time_elapsed/60,
                time_remaining, time_remaining/60

            ))

            rep_pos += 1

    for key in result_dict.keys():

        filename = os.path.join(support_folder, base_filename.format(key=key, a=a_curr, i=i_curr))
        np.save(filename, result_dict[key])

    time_elapsed = time.time() - sub_start

    print("{} - Done! -> {:.0f}s ({:.0f}min) elapsed".format(current_name, time_elapsed, time_elapsed / 60))


if __name__ == '__main__':

    print("Calculation Started!!")

    with concurrent.futures.ProcessPoolExecutor() as executor:

        future_list = list()

        for a in range(len(p_rel_list)):

            for i in range(len(t_rel_list)):

                time.sleep(0.5)
                future_list.append(executor.submit(evaluate_points, a, i))

        elapsed_time_list = list()
        for f in concurrent.futures.as_completed(future_list):
            elapsed_time_list.append(f.result())

    print("Calculation Completed!!")
    print("Saving Results ...")
    time.sleep(2)
    base_shape = np.array(grad_nd_nd.shape)
    res_shape = np.append([len(p_rel_list), len(t_rel_list)], base_shape)
    ovr_res_dict = init_result_dict(res_shape)

    for key in ovr_res_dict.keys():

        for a in range(len(p_rel_list)):

            for i in range(len(t_rel_list)):

                filename = os.path.join(support_folder, base_filename.format(

                    key=key, a=a, i=i

                ))

                if os.path.exists(filename):

                    ovr_res_dict[key][a, i, :, :] = np.load(filename)
                    # os.remove(filename)

        filename = os.path.join(result_folder, "{}.npy".format(key))
        np.save(filename, ovr_res_dict[key])

    filename = os.path.join(result_folder, "elapsed_times.npy")
    np.save(filename, np.array(elapsed_time_list))

    filename = os.path.join(result_folder, "dz_nd.npy")
    np.save(filename, dz_nd)

    filename = os.path.join(result_folder, "grad_nd_nd.npy")
    np.save(filename, grad_nd_nd)

    print("Calculation Completed!!!")

