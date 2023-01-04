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

real T={max_time};
int nT={time_steps};
real deltat=T/nT;

real gradRock={grad_rock}, DT={DT};"""

__PLOT_AND_EXPORT_RESULTS = """
    if (counter%plotevery == 0){{

        plot(uInt, value=true, grey=true, fill=true);   

    }}

    if (counter%saveevery == 0){{

        sum = 0;
        fx = dx(uInt);
        fy = dy(uInt);

        for (int i = 0; i < 1000; i++){{

            theta = pi * i / 1000;
            x = r*sin(theta);
            y = r*cos(theta);
            sum = sum + (fx * sin(theta) + fy * cos(theta));

        }}
        
        ff << t << "; " << sum/1000 * 2*pi*r << endl;   
        
    }}
"""

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

// field initializer
func u0 = -gradRock*y;

// The finite element space defined over ThInt is called here VhInt
fespace VhInt(ThInt, P1);
VhInt uInt, uIntOld, uIntPrev=u0, vInt;	// Define uInt and vInt as piecewise-P1 continuous functions

// The finite element space defined over ThInt is called here VhInt
fespace VhExt(ThExt, P1);
VhExt uExt, vExt, uExtPrev=u0, g=gradRock;	// Define uInt and vInt as piecewise-P1 continuous functions

// DEFINE INTERNAL PROBLEM
problem thermalInt(uInt, vInt, solver=LU)

	= int2d(ThInt)(    

        uInt * vInt
		+ deltat * (

		    dx(uInt)*dx(vInt)	
		    + dy(uInt)*dy(vInt)

		)

	) - int2d(ThInt)(

        uIntPrev * vInt

	) + on(Ctube, uInt=-DT) + on(CMeshInt, uInt=uExtPrev);
	
problem thermalIntNOFLUX(uInt, vInt, solver=LU)

	= int2d(ThInt)(    

        uInt * vInt
		+ deltat * (

		    dx(uInt)*dx(vInt)	
		    + dy(uInt)*dy(vInt)

		)

	) - int2d(ThInt)(

        uIntPrev * vInt

	) + on(CMeshInt, uInt=uExtPrev);

// DEFINE EXTERNAL PROBLEM
problem thermalExt(uExt, vExt, solver=LU)

	= int2d(ThExt)(    // The bilinear part

        uExt*vExt
		+ deltat * (

		    dx(uExt)*dx(vExt)	
		    + dy(uExt)*dy(vExt)

		)

	) - int2d(ThExt)(

        uExtPrev * vExt

	) - int1d(ThExt, CMeshExt)(

		-(dx(uIntPrev)*x/(ratioR*r)+ dy(uIntPrev)*y/(ratioR*r))*vExt*deltat

	) - int1d(ThExt, Cbottom)(

        g*vExt*deltat

    )
	+ on(Ctop, uExt=-gradRock*H);

real t=0;
int counter=0;

int plotevery = nT / 50;
int saveevery = nT / 200;
real theta, sum=0;
VhInt fx, fy;

ofstream ff("thermic.dat");

while(t < T){{
    
    if (t%(T/3) < T/6){{
    
	    thermalInt;
	
	}} else {{
	
	    thermalIntNOFLUX;
	    
	}}
	
	thermalExt;

    {plot_export_result}

    t = t + deltat;
    counter = counter + 1;

    uIntPrev=uInt;
    uExtPrev=uExt;

}}

"""

def write_edp_file(

    n_points=10, n_points_circle=30,
    ratio_L=2000, ratio_H=1000,
    DT=1, grad_rock=0.01,
    max_time_year=5, time_steps=10000

):

    param_init = __MAIN_PARAMETER_INITIALIZATION.format(

        n_points=n_points, n_points_circle=n_points_circle,
        ratio_L=ratio_L, ratio_H=ratio_H,
        DT=DT, grad_rock=grad_rock,

        max_time=max_time_year*365*24*3600*1E-6, time_steps=time_steps

    )

    with open(EDP_FILE_PATH, "w") as edp_file:

        edp_file.write(

            __DOUBLE_MESH_CALCULATIONS.format(

                main_parameter_initialization=param_init,
                plot_export_result=__PLOT_AND_EXPORT_RESULTS

            )

        )


if __name__ == "__main__":

    write_edp_file()