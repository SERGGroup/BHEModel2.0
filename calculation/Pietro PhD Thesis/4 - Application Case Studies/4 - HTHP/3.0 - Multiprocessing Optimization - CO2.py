from main_code.power_plants.HTHP.subclasses.CO2_heat_pump_py import CO2HeatPumpThermo
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.constants import CALCULATION_FOLDER, os
from scipy.optimize import minimize
import concurrent.futures
import numpy as np
import math
import time
import gc


result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "4 - Application Case Studies", "4 - HTHP",
    "res"

)

support_folder = os.path.join(result_folder, "support")

base_filename = "res_{i}_{j}.npy"
reporting_points = [0.05, 0.25, 0.5, 0.75, 1]

p_steam = 2
t_in_wells = [10, 50]
omegas = [0.5, 1, 2]
n_temp = len(t_in_wells)

n_points = 75
depths_list = np.logspace(np.log10(250), np.log10(5000), n_points)
grads_list = np.logspace(np.log10(25), 2, n_points)
depths, grads = np.meshgrid(depths_list, grads_list, indexing="ij")

fluid = "Carbon Dioxide"
bhe_in = PlantThermoPoint([fluid], [1])
bhe_in.set_variable("T", bhe_in.RPHandler.TC)
bhe_in.set_variable("P", bhe_in.RPHandler.PC)
rho_crit = bhe_in.get_variable("rho")


def update_values(depth, t_rock, t_amb, p_in, t_sc_perc):

    new_hthp = CO2HeatPumpThermo(

        P_steam=p_steam,
        BHE_depth=depth, T_rock=t_rock,
        T_in_BHE=t_amb, P_in_BHE=p_in,
        T_ambient=t_amb - 10

    )

    new_hthp.T_sc_perc = t_sc_perc
    new_hthp.calculate(calculate_thermo_only=True)

    return new_hthp


def evaluate_part(i, j):

    rep_pos = 0
    curr_value = i * n_points + j + 1
    max_digits = math.ceil(np.log10(n_points*n_temp))
    current_name = "{{:{}d}} of {{}}".format(max_digits).format(curr_value, n_points*n_temp)

    print("{} - Started".format(current_name))

    sub_start = time.time()

    sub_results = np.empty([4 * len(omegas), n_points])
    sub_results[:] = np.nan

    bhe_in.set_variable("T", t_in_wells[i])
    bhe_in.set_variable("rho", rho_crit)
    p_ref = bhe_in.get_variable("P")

    for k in range(n_points):

        t_rocks = t_in_wells[i] + grads[j, k] * depths[j, k] / 1000

        def opt_func(x, omega=(1,)):

            p_curr = p_ref * x[0]
            base_hthp = update_values(

                depth=depths[j, k],
                t_rock=t_rocks,
                t_amb=t_in_wells[i],
                p_in=p_curr,
                t_sc_perc=x[1]

            )

            return 1 / (20 * base_hthp.m_dot_ratio) + omega * 3 / base_hthp.COP

        x0 = [2, 0.5]
        for n in range(len(omegas)):

            res = minimize(opt_func, np.array(x0), args=[n], bounds=[(None, None), (0, 1)])
            opt_res = -res.fun
            opt_res_x = res.x[0]

            res_test = -opt_func([0.9999, res.x[1]], omega=(n,))
            if res_test > opt_res:
                opt_res_x = [0.9999, res.x[1]]
                opt_res = res_test

            res_test = -opt_func([1.00001, res.x[1]], omega=(n,))
            if res_test > opt_res:
                opt_res_x = [1.00001, res.x[1]]

            opt_hthp = update_values(

                depth= depths[j, k],
                t_rock= t_rocks,
                t_amb= t_in_wells[i],
                p_in= p_ref * opt_res_x[0],
                t_sc_perc= opt_res_x[1]

            )

            sub_results[n, k] = opt_res_x[0]
            sub_results[n + 3, k] = opt_res_x[1]
            sub_results[n + 6, k] = opt_hthp.COP
            sub_results[n + 9, k] = opt_hthp.m_dot_ratio

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

            elapsed_time_list.append(f.result())

    print("Calculation Completed!!")
    print("Saving Results ...")
    time.sleep(2)

    res_shape = np.append([4 * len(omegas), n_temp], grads.shape)
    results = np.empty(res_shape)
    results[:] = np.nan

    for i in range(n_temp):

        for j in range(n_points):

            filename = os.path.join(support_folder, base_filename.format(i=i, j=j))

            if os.path.exists(filename):

                results[:, i, j, :] = np.load(filename)
                # os.remove(filename)

    filename = os.path.join(result_folder, "opt_results.npy")
    np.save(filename, results)

    filename = os.path.join(result_folder, "depths_mesh.npy")
    np.save(filename, depths)

    filename = os.path.join(result_folder, "grads_mesh.npy")
    np.save(filename, grads)

    filename = os.path.join(result_folder, "t_ambs.npy")
    np.save(filename, t_in_wells)

    filename = os.path.join(result_folder, "elapsed_times.npy")
    np.save(filename, np.array(elapsed_time_list))

    print("Calculation Completed!!!")
