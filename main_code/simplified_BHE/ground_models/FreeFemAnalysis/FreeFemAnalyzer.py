from BHEOther.BHESimplifiedAnalysis.new_analysis.code.simplified_BHE.ground_models.FreeFemAnalysis.edp_costants import (

    EDP_FILE_PATH, FREE_FEM_FOLDER, FREE_FEM_EXEC,
    write_edp_file

)
import matplotlib.pyplot as plt
import os, subprocess


def retrieve_result():

    output_file_path = os.path.join(FREE_FEM_FOLDER, "mean_gradient.txt")

    if os.path.isfile(output_file_path):

        with open(output_file_path, "r") as edp_file:

            line = edp_file.readline()

        os.remove(output_file_path)
        return float(line)

    else:

        return None


def calculate(

        temperature_bc=False, perfrom_double_mesh_calculation=True,
        overall_mesh=True, mesh_refinement=False,
        n_points=7, n_points_circle=20,
        ratio_L=2000, ratio_H=1000,
        DT=40., grad_rock=0.1,
        grad_tube_ratio=10

):

    write_edp_file(

        temperature_bc=temperature_bc, perfrom_double_mesh_calculation=perfrom_double_mesh_calculation,
        overall_mesh=overall_mesh, mesh_refinement=mesh_refinement,
        n_points=n_points, n_points_circle=n_points_circle,
        ratio_L=ratio_L, ratio_H=ratio_H,
        DT=DT, grad_rock=grad_rock,
        grad_tube_ratio=grad_tube_ratio

    )

    subprocess.call(

        '"{}" "{}"'.format(FREE_FEM_EXEC, EDP_FILE_PATH),
        cwd=FREE_FEM_FOLDER, shell=True,
        stdout=open(os.devnull, 'w'),
        stderr=subprocess.STDOUT

    )

    return retrieve_result()


if __name__ == "__main__":

    fig, ax = plt.subplots()

    n_points = 10
    for grad_rock in [0.01]:

        x = list()
        y = list()

        for i in range(100):

            ratio_L = 10 ** 4
            ratio_Q = 1 + i

            x.append(ratio_Q)
            y.append(calculate(

                ratio_L=ratio_L,
                ratio_H=int(ratio_L/2),
                temperature_bc=False,
                perfrom_double_mesh_calculation=True,
                n_points=n_points, n_points_circle=n_points*10,
                grad_tube_ratio=ratio_Q, grad_rock=grad_rock

            ))

            print("{}, {} -> {}".format(grad_rock, ratio_Q, y[-1]))

        ax.plot(x, y, label=str(grad_rock))

    # plt.xscale("log")
    ax.legend()
    plt.show()