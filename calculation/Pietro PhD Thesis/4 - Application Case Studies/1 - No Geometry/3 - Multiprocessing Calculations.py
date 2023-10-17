from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.power_plants.base_surface_power_plant import BaseSurfacePlant
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.constants import CALCULATION_FOLDER, os
from scipy.optimize import minimize
import multiprocessing
import numpy as np
import time
import gc

result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "4 - Application Case Studies", "1 - No Geometry",
    "res", "support"

)

base_filename = "res_{i}_{j}.npy"

n_temp = 3
t_ambs = np.linspace(10, 50, n_temp)

n_points = 25
depths_list = np.logspace(2, 4, n_points)
grads_list = np.logspace(1, 2, n_points)
depths, grads = np.meshgrid(depths_list, grads_list, indexing="ij")

fluid = "Carbon Dioxide"
bhe_in = PlantThermoPoint([fluid], [1])
bhe_in.set_variable("T", bhe_in.RPHandler.TC)
bhe_in.set_variable("P", bhe_in.RPHandler.PC)
rho_crit = bhe_in.get_variable("rho")

def update_values(curr_bhe, t_amb, p_in):

    curr_bhe.input_point.set_variable("T", t_amb)
    curr_bhe.input_point.set_variable("P", p_in)

    new_base_surface = BaseSurfacePlant(curr_bhe)
    new_base_surface.prepare_calculation()
    new_base_surface.calculate()

    return new_base_surface

def evaluate_part(i, j):

    print(f'{i}-{j} Started')
    sub_start = time.time()

    sub_results = np.empty([2 * 3, n_points])
    sub_results[:] = np.nan

    bhe_in.set_variable("T", t_ambs[i])
    bhe_in.set_variable("rho", rho_crit)
    p_ref = bhe_in.get_variable("P")

    for k in range(n_points):

        t_rocks = t_ambs[i] + grads[j, k] * depths[j, k] / 1000
        bhe = SimplifiedBHE(bhe_in, dz_well=depths[j, k], t_rocks=t_rocks)

        def opt_func(x, opt_value=(0,)):

            p_curr = p_ref * x[0]
            base_surface = update_values(curr_bhe=bhe, t_amb=t_ambs[i], p_in=p_curr)

            if opt_value[0] == 0:

                return -base_surface.w_dot

            elif opt_value[0] == 1:

                return -base_surface.ex_dot

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

    filename = os.path.join(result_folder, base_filename.format(i=i, j=j))
    np.save(filename, sub_results)

    gc.collect()
    print(f'{i}-{j} Done - {time.time() - sub_start}s')


if __name__ == '__main__':

    print("Calculation Started!!")

    processes = list()
    for i in range(n_temp):

        for j in range(n_points):

            p = multiprocessing.Process(target=evaluate_part, args=[i, j])
            p.start()

            processes.append(p)

    for process in processes:

        process.join()

    print("Calculation Completed!!")
    print("Saving Results ...")
    time.sleep(2)

    res_shape = np.append([2 * 3, n_temp], grads.shape)
    results = np.empty(res_shape)
    results[:] = np.nan

    for i in range(n_temp):

        for j in range(n_points):

            filename = os.path.join(result_folder, base_filename.format(i=i, j=j))

            if os.path.exists(filename):

                results[:, i, j, :] = np.load(filename)
                os.remove(filename)

    print(results)
    filename = os.path.join(result_folder, "Optimization Result")
    np.save(filename, results)
