# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.well_model.simplified_well.simplified_well import SimplifiedBHE
from main_code.support.abstract_plant_thermo_point import PlantThermoPoint
from main_code.support.other.excel_exporter import write_excel_sheet
from main_code.constants import CALCULATION_FOLDER, os
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   DEFINE WORKING FOLDERS                ------------------------------------->
result_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "2 - Theoretical Analysis", "res"

)


# %%-------------------------------------   RESULT DF FUNCTIONS DEFINITION      -------------------------------------> #
def init_result_df():

    return {

        'dz': {"unit": ["m"], "values": [list()]},
        'geo_grad': {"unit": ["°C/m"], "values": [list()]},
        'p_rel': {"unit": ["-"], "values": [list()]},
        't_sat_mult': {"unit": ["-"], "values": [list()]},

        'beta': {"unit": ["-"], "values": [list()]},
        'tau': {"unit": ["-"], "values": [list()]},

        'dh_exp [kJ/kg]': {"unit": ["Ovrll", "Exp", "Cool"], "values": [list(), list(), list()]},
        'dex_exp [kJ/kg]': {"unit": ["Ovrll", "Exp", "Cool"], "values": [list(), list(), list()]},

        'dh_cool [kJ/kg]': {"unit": ["Ovrll", "Exp", "Cool"], "values": [list(), list(), list()]},
        'dex_cool [kJ/kg]': {"unit": ["Ovrll", "Exp", "Cool"], "values": [list(), list(), list()]},

    }

def append_results(res_dict:dict, dz_curr, grad_curr, p_rel_curr, t_sat_mult_curr, bhe_res: SimplifiedBHE):

    if not (bhe_res.input_point.get_variable("rho") < 0 or bhe_res.output_point.get_variable("h") < -1e6):

        res_dict["dz"]["values"][0].append(dz_curr)
        res_dict["geo_grad"]["values"][0].append(grad_curr)
        res_dict["p_rel"]["values"][0].append(p_rel_curr)
        res_dict["t_sat_mult"]["values"][0].append(t_sat_mult_curr)

        beta = bhe_res.output_point.get_variable("P") / bhe_res.input_point.get_variable("P")

        if bhe_res.input_point.get_unit("T") == "K":
            tau = bhe_res.t_rocks / bhe_res.input_point.get_variable("T")

        else:
            tau = (bhe_res.t_rocks + 273.15) / (bhe_res.input_point.get_variable("T") + 273.15)

        res_dict["beta"]["values"][0].append(beta)
        res_dict["tau"]["values"][0].append(tau)

        output_point = bhe_res.output_point
        input_point = bhe_res.input_point

        exp_first_point = bhe_res.output_point.duplicate()
        exp_first_point.set_variable("P", bhe_res.input_point.get_variable("P"))
        exp_first_point.set_variable("s", bhe_res.output_point.get_variable("s"))

        res_dict['dh_exp [kJ/kg]']["values"][0].append(output_point.dvar(input_point, "h"))
        res_dict['dh_exp [kJ/kg]']["values"][1].append(output_point.dvar(exp_first_point, "h"))
        res_dict['dh_exp [kJ/kg]']["values"][2].append(exp_first_point.dvar(input_point, "h"))

        t_ref = input_point.get_variable("T")
        res_dict['dex_exp [kJ/kg]']["values"][0].append(output_point.dex(input_point, t_ref))
        res_dict['dex_exp [kJ/kg]']["values"][1].append(output_point.dex(exp_first_point, t_ref))
        res_dict['dex_exp [kJ/kg]']["values"][2].append(exp_first_point.dex(input_point, t_ref))

        cool_first_s_point = bhe_res.output_point.duplicate()
        cool_first_s_point.set_variable("P", bhe_res.output_point.get_variable("P"))
        cool_first_s_point.set_variable("s", bhe_res.input_point.get_variable("s"))
        res_dict['dh_cool [kJ/kg]']["values"][0].append(output_point.dvar(input_point, "h"))
        res_dict['dh_cool [kJ/kg]']["values"][1].append(cool_first_s_point.dvar(input_point, "h"))
        res_dict['dh_cool [kJ/kg]']["values"][2].append(output_point.dvar(cool_first_s_point, "h"))

        res_dict['dex_cool [kJ/kg]']["values"][0].append(output_point.dex(input_point, t_ref))
        res_dict['dex_cool [kJ/kg]']["values"][1].append(cool_first_s_point.dex(input_point, t_ref))
        res_dict['dex_cool [kJ/kg]']["values"][2].append(output_point.dex(cool_first_s_point, t_ref))

def extract_sub_dict(res_dict:dict, key, value, sub_list=0):

    indices_list = list(

        filter(

            lambda x: res_dict[key]["values"][sub_list][x] == value,
            range(len(res_dict[key]["values"][sub_list]))

        )

    )

    new_dict = init_result_df()

    for index in indices_list:

        for key in res_dict.keys():

            for i in range(len(new_dict[key]["values"])):

                new_dict[key]["values"][i].append(res_dict[key]["values"][i][index])

    return new_dict


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
input_point = PlantThermoPoint(["CarbonDioxide"], [1])
p_crit = input_point.RPHandler.PC

dz_list = np.logspace(2, 4, 5)
grad_list = np.linspace(10, 100, 5) / 1e3
p_rels = [10 ** -4, 10 ** -2, 10 ** 0]
t_sat_mults = [0.8, 0.9, 1.1, 2, 4]

res_dict = init_result_df()
pbar = tqdm(desc="Calculating Points", total=len(dz_list)*len(grad_list)*len(p_rels)*len(t_sat_mults))

for dz in dz_list:

    for grad in grad_list:

        for p_rel in p_rels:

            p_in = p_rel * p_crit
            input_point.set_variable("P", p_in)
            input_point.set_variable("Q", 0)
            t_sat = input_point.get_variable("T")

            for t_sat_mult in t_sat_mults:

                t_in = t_sat_mult * t_sat
                input_point.set_variable("P", p_in)
                input_point.set_variable("T", t_in)

                bhe_CO2 = SimplifiedBHE(

                    input_thermo_point=input_point,
                    dz_well=dz, t_rocks=t_in + grad * dz

                )

                bhe_CO2.update()

                append_results(

                    res_dict,
                    dz_curr=dz, grad_curr=grad,
                    p_rel_curr=p_rel, t_sat_mult_curr=t_sat_mult,
                    bhe_res=bhe_CO2

                )

                pbar.update(1)

pbar.close()

# %%-------------------------------------   SAVE RESULTS                        -------------------------------------> #
file_path = os.path.join(result_folder, "1 - base comparison.xlsx")
write_excel_sheet(excel_path=file_path, sheet_name="results", data_frame=res_dict, overwrite="hard")


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #

result = extract_sub_dict(res_dict, "p_rel", 10 ** -2)

fig, axs = plt.subplots(1, len(t_sat_mults), dpi=300)
fig.set_size_inches(5 * len(t_sat_mults), 4)

for i in range(len(t_sat_mults)):

    ax = axs[i]
    ax.set_title("t_sat_mult = {}".format(t_sat_mults[i]))

    sub_res = extract_sub_dict(result, "t_sat_mult", t_sat_mults[i])

    for grad in grad_list:

        sub_sub_res = extract_sub_dict(sub_res, "geo_grad", grad)
        ax.plot(sub_sub_res["dz"]["values"][0], sub_sub_res["tau"]["values"][0], label="geo_grad={}".format(grad))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()

plt.tight_layout()
plt.show()


# %%-------------------------------------   PLOT NON-DIMENSIONAL RESULT         -------------------------------------> #
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
styles = ["+", "x", "1", "2", "3", "4"]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

for i in range(len(t_sat_mults)):

    sub_res = extract_sub_dict(res_dict, "t_sat_mult", t_sat_mults[i])

    for j in range(len(p_rels)):

        sub_sub_res = extract_sub_dict(sub_res, "p_rel", p_rels[j])
        ax.scatter(

            np.array(sub_sub_res["tau"]["values"][0]) - 1,
            np.array(sub_sub_res["tau"]["values"][0]) - 1,
            sub_sub_res["dex_exp [kJ/kg]"]["values"][0],
            marker=styles[i], color=colors[j]

        )

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_zscale("log")

ax.set_xlim((10**-4, 10**1))
plt.tight_layout()
plt.show()