from constants import FREE_FEM_FOLDER
import os

EDP_FILE_PATH = os.path.join(FREE_FEM_FOLDER, "Temperature Distribution.edp")
FREE_FEM_EXEC = "C:\Program Files (x86)\FreeFem++\FreeFem++.exe"

# <------------------------------------------------------------------------->
# <------------------------------------------------------------------------->
# <-----------------------   OVERALL PARAMETERS  --------------------------->
# <------------------------------------------------------------------------->
# <------------------------------------------------------------------------->

__MAIN_PARAMETER_INITIALIZATION = """
int pointsH={n_points}, pointsCircle={n_points_circle};
real ratioL={ratio_L}, ratioH={ratio_H};
real r=1., L=r*ratioL, H=r*ratioH;
int pointsL = pointsH*L/H;

real gradRock={grad_rock},"""

__EXPORT_RESULTS = """
int K = 1000; 
real t, sum=0;

for (int i = 0; i < K; i++){{
	
	t = pi * i / K;
    x = r*sin(t);
	y = r*cos(t);
	sum = sum + {result_element};

}}

real tubeFlux = sum / K * 2 * pi * r;
ofstream file("mean_gradient.txt");
file.precision(10);
file << tubeFlux;
"""

__TUBE_TEMP_BC = {

    "param_init": __MAIN_PARAMETER_INITIALIZATION + " DT={DT};",
    "bc_definition": """ + on(Ctube, uInt=-DT)""",
    "result_export": __EXPORT_RESULTS.format(result_element="fx * sin(t) + fy * cos(t)")
}

__TUBE_GRAD_BC = {

    "param_init": __MAIN_PARAMETER_INITIALIZATION + """ gradTubeRatio={grad_tube_ratio};
real gradTube=gradRock*gradTubeRatio;
    """,
    "bc_definition": """- int1d(ThInt, Ctube)(
    
		-gradTube*vInt
	
	)""",
    "result_export": __EXPORT_RESULTS.format(result_element=" u")
}

def write_edp_file(

    temperature_bc=True, perfrom_double_mesh_calculation=True,
    overall_mesh=True, mesh_refinement=False,
    n_points=7, n_points_circle=20,
    ratio_L=2000, ratio_H=1000,
    DT=20., grad_rock=0.01,
    grad_tube_ratio=1

):

    if temperature_bc:

        param_init = __TUBE_TEMP_BC["param_init"].format(

            n_points=n_points, n_points_circle=n_points_circle,
            ratio_L=ratio_L, ratio_H=ratio_H,
            DT=DT, grad_rock=grad_rock,

        )
        export_result = __TUBE_TEMP_BC["result_export"]
        bc_definition = __TUBE_TEMP_BC["bc_definition"]

    else:

        param_init = __TUBE_GRAD_BC["param_init"].format(

            n_points=n_points, n_points_circle=n_points_circle,
            ratio_L=ratio_L, ratio_H=ratio_H,
            grad_tube_ratio=grad_tube_ratio, grad_rock=grad_rock,

        )
        export_result = __TUBE_GRAD_BC["result_export"]
        bc_definition = __TUBE_GRAD_BC["bc_definition"]

    if perfrom_double_mesh_calculation:

        __write_double_calculation_edp(

            param_init, export_result, bc_definition

        )

    else:

        __write_simple_calculation_edp(

            param_init, export_result=export_result, bc_definition=bc_definition,
            overall_mesh=overall_mesh, mesh_refinement=mesh_refinement

        )

# <------------------------------------------------------------------------->
# <------------------------------------------------------------------------->
# <-------------------   SIMPLE MESH CALCULATIONS   ------------------------>
# <------------------------------------------------------------------------->
# <------------------------------------------------------------------------->

__SIMPLE_CALCULATION_STRING = """
// The finite element space defined over Th is called here Vh
fespace Vh(Th, P1);
Vh u, v, g=gradRock;	// Define u and v as piecewise-P1 continuous functions

// Define the PDE
solve thermal(u, v, solver=LU)

    = int2d(Th)(// The bilinear part

        dx(u)*dx(v)	
        + dy(u)*dy(v)

    ) - int1d(Th, Cbottom)(

        g*v

    ) {tube_boundary_condition} 
    + on(Ctop, u=-gradRock*H);	// The Dirichlet boundary condition

// Evaluate the mean gradient on the internal circle
Vh fx, fy;
fx = dx(u);
fy = dy(u);

{export_results_script}
"""

__EDP_SCRIPT_STRING = """
include "ffmatlib.idp"

// Parameters
{main_parameter_initialization}

// Define mesh boundary
int Cside=99, Ctop=98, Cbottom=97, Ctube=96;
real refstart=0.10*(H - r) + r;

border C1(t=2*L, 0){{x=t; y=H; label=Ctop;}}
border C21(t=H, refstart){{x=0; y=t; label=Cside;}}
border C22(t=refstart, r){{x=0; y=t; label=Cside;}}
border C23(t=0, pi){{x=r*sin(t); y=r*cos(t); label=Ctube;}}
border C24(t=-r, -refstart){{x=0; y=t; label=Cside;}}
border C25(t=-refstart, -H){{x=0; y=t; label=Cside;}}
border C3(t=0, 2*L){{x=t; y=-H; label=Cbottom;}}
border C4(t=-H, H){{x=2*L; y=t; label=Cside;}}

// The triangulated domain Th is on the left side of its boundary
mesh Th = buildmesh(C1(pointsL)+C21(pointsH/2)+C22(pointsCircle)+C23(pointsCircle)+C24(pointsCircle)+C25(pointsH/2)+C3(pointsL)+C4(pointsH));
{optional_mesh_refinement}
{calculation_string}

"""

__EDP_SCRIPT_OVERALL_STRING = """
include "ffmatlib.idp"

// Parameters
{main_parameter_initialization}

// Define mesh boundary
int Cside=99, Ctop=98, Cbottom=97, Ctube=96;

border C1(t=L, -L){{x=t; y=H; label=Ctop;}}
border C2(t=H, -H){{x=-L; y=t; label=Cside;}}
border C3(t=-L, L){{x=t; y=-H; label=Cbottom;}}
border C4(t=-H, H){{x=L; y=t; label=Cside;}}
border C0(t=0, 2*pi){{x=r*sin(t); y=r*cos(t); label=Ctube;}}

// The triangulated domain Th is on the left side of its boundary
mesh Th = buildmesh(C0(pointsCircle)+C1(pointsL)+C2(pointsH)+C3(pointsL)+C4(pointsH));

{optional_mesh_refinement}
{calculation_string}

"""

__OPTIONAL_MESH_REFINEMENT = """

int refIndex=5;
Th = splitmesh(Th, 2 + int(refIndex * r / sqrt(x*x + y*y)));

"""

def __write_simple_calculation_edp(

    param_init, export_result, bc_definition,
    overall_mesh, mesh_refinement

):
    if mesh_refinement:

        refinement_str = __OPTIONAL_MESH_REFINEMENT

    else:

        refinement_str = ""

    if overall_mesh:

        base_str = __EDP_SCRIPT_OVERALL_STRING

    else:

        base_str = __EDP_SCRIPT_STRING

    with open(EDP_FILE_PATH, "w") as edp_file:

        edp_file.write(

            base_str.format(

                main_parameter_initialization=param_init,
                optional_mesh_refinement=refinement_str,
                calculation_string=__SIMPLE_CALCULATION_STRING.format(

                    export_results_script=export_result,
                    tube_boundary_condition=bc_definition

                )

            )

        )

# <------------------------------------------------------------------------->
# <------------------------------------------------------------------------->
# <-------------------   DOUBLE MESH CALCULATIONS   ------------------------>
# <------------------------------------------------------------------------->
# <------------------------------------------------------------------------->

__DOUBLE_MESH_CALCULATIONS = """

include "ffmatlib.idp"

// Parameters
{main_parameter_initialization}

int Cside=99, Ctop=98, Cbottom=97,CMeshExt=96, CMeshInt=95, Ctube=94;
real ratioR=100;

// external mesh boundary
border C1(t=L, -L){{x=t; y=H; label=Ctop;}}
border C2(t=H, -H){{x=-L; y=t; label=Cside;}}
border C3(t=-L, L){{x=t; y=-H; label=Cbottom;}}
border C4(t=-H, H){{x=L; y=t; label=Cside;}}
border C10(t=0, 2*pi){{x=ratioR*r*sin(t); y=ratioR*r*cos(t); label=CMeshExt;}}

// internal mesh boundary
border C0(t=0, 2*pi){{x=r*sin(t); y=r*cos(t); label=Ctube;}}
border C01(t=2*pi, 0){{x=ratioR*r*sin(t); y=ratioR*r*cos(t); label=CMeshInt;}}

// The triangulated domain Th is on the left side of its boundary
mesh ThInt = buildmesh(C0(pointsCircle)+C01(pointsCircle));
mesh ThExt = buildmesh(C10(pointsCircle)+C1(pointsL)+C2(pointsH)+C3(pointsL)+C4(pointsH));

// The finite element space defined over ThInt is called here VhInt
fespace VhInt(ThInt, P1);
VhInt uInt, uIntOld, vInt;	// Define uInt and vInt as piecewise-P1 continuous functions

// The finite element space defined over ThInt is called here VhInt
fespace VhExt(ThExt, P1);
VhExt uExt, vExt, g=gradRock;	// Define uInt and vInt as piecewise-P1 continuous functions

// DEFINE INTERNAL PROBLEM
problem thermalInt(uInt, vInt, solver=LU)

	= int2d(ThInt)(    // The bilinear part

		dx(uInt)*dx(vInt)	
		+ dy(uInt)*dy(vInt)

	) {tube_boundary_condition} + on(CMeshInt, uInt=uExt);	// The Dirichlet boundary condition

// DEFINE EXTERNAL PROBLEM
problem thermalExt(uExt, vExt, solver=LU)

	= int2d(ThExt)(    // The bilinear part

		dx(uExt)*dx(vExt)	
		+ dy(uExt)*dy(vExt)

	) - int1d(ThExt, CMeshExt)(

		-(dx(uInt)*x/(ratioR*r)+ dy(uInt)*y/(ratioR*r))*vExt

	) - int1d(ThExt, Cbottom)(

        g*vExt

    )
	+ on(Ctop, uExt=-gradRock*H);	// The Dirichlet boundary condition

func real error() {{

	real t, fnet, sum=0, sumExt=0;

	for (int i = 0; i < 500; i++){{

		x = ratioR*r*sin(t);
		y = ratioR*r*cos(t);
		sum = sum + abs(uInt - uExt);
		sumExt = sumExt + abs(uExt);
	}}

	return(sum / sumExt);

}}

int n=0;
real alpha={alpha};
while(n < {max_iterations}){{

	uIntOld = uInt;

	thermalInt;
	uInt = alpha*uIntOld + (1-alpha)*uInt;
	thermalExt;

	if (error() < {toll}) break;

}}

// Evaluate the mean gradient on the internal circle
VhInt fx, fy, u;
fx = dx(uInt);
fy = dy(uInt);
u = uInt;

{export_results_script}

"""

def __write_double_calculation_edp(

    param_init, export_result, bc_definition,
    alpha=0.5, max_iterations=50,
    toll = 1E-5

):

    with open(EDP_FILE_PATH, "w") as edp_file:

        edp_file.write(

            __DOUBLE_MESH_CALCULATIONS.format(

                alpha=alpha, max_iterations=max_iterations, toll=toll,
                main_parameter_initialization=param_init,
                export_results_script=export_result,
                tube_boundary_condition=bc_definition

            )

        )