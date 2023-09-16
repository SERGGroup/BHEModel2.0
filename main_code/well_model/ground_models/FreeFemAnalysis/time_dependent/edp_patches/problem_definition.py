TIME_ITERATION_PARAM = """
real Tstart={min_time};
real Tend={max_time};
int nT={time_steps};
{delta_t_definition}
"""

STD_DELTAT = {

    "Definition":"""
real deltat = Tend / nt;
""",

    "Implementation": """
    t = deltat * i;
"""
}

LOG_DELTAT = {

    "Definition": """
real logdeltat = (log10(Tend)-log10(Tstart)) / nT;
real deltat = Tstart;
""",

    "Implementation": """
    t = 10^(logdeltat * i + log10(Tstart));
    deltat = t - told;
"""
}

PROBLEM_DEF_PARAM = """
// Problem Definition Parameters
real gradRock={grad_rock}, DT={DT};
real vxUnd={vx}, vyUnd={vy};
real Pe={Pe};
"""

SPEED_ADIMENSIONALIZATION = """
// Speed Adimensionalization
real modvUnd = sqrt(vxUnd ^ 2 + vyUnd^2);
vxUnd = vxUnd / modvUnd;
vyUnd = vyUnd / modvUnd;

// Speed Field Definition 
// (Potential Velocity Field)
real cosAlphaV = vxUnd;
real sinAlphaV = vyUnd;
func cosThetaV = (x * vxUnd + y * vyUnd) / sqrt(x^2 + y^2);
func sinThetaV = sqrt(1 - cosThetaV^2);

func rInv = 1 / sqrt(x^2 + y^2);
func potFieldFunc = (1/rInv + rInv) * cosThetaV;
"""
EMPTY_FIELD_INIT = """
func potFieldFunc = 0;
"""

PROBLEM_DEFINITION = """
// The finite element space defined over Th is called here Vh
fespace Vh(Th, P1);
func u0 = -gradRock*y;
Vh u, v, uPrev=u0, f=gradRock;	// Define u and v as piecewise-P1 continuous functions
Vh potField=potFieldFunc;

// The finite element space defined over ThPlot for visualization only
fespace Vplot(ThPlot, P1);
Vplot uplot = uPrev;

// Define the PDE
solve thermal(u, v, solver=LU)

	= int2d(Th)(    // The bilinear part
		u * v
		+ deltat * (
		
            dx(u)*dx(v)	
            + dy(u)*dy(v) +
            Pe * (dx(potField) * dx(u) + dy(potField) * dy(u)) * v
        
        )
		
	) - int2d(Th)(
	
        uPrev * v

	) - int1d(Th, Cbottom)(
    
		f*v*deltat
	
	)
	+ on(Ctop, u=-gradRock*H)
	+ on(Ctube, u=-DT);	   // The Dirichlet boundary condition

"""

SOLUTION_AND_SAVING = """
real told=0, t=0;
int nSections={n_sections_save};

{plot_every_initialization}
int saveevery = nT / {n_save};
real theta, sum=0;
Vh fx, fy;

{open_save_stream}

for (int i = 0; i <= nT; ++i){{

    told = t;
    {delta_t_implementation}

	thermal;

    {display_plot}
    {save_output}

    uPrev=u;

}}
"""

DISPLAY_PLOTS = """
    if (i%plotevery == 0){{

        uplot = u;
        plot(uplot, value=true, grey=true, fill=true);   

    }}
"""

PLOT_EVERY_INIT = """
int plotevery = nT / {n_plot};
"""

SAVE_OUTPUT = """
    if (i%saveevery == 0){{

        sum = 0;
        fx = dx(u);
        fy = dy(u);

        for (int i = 0; i < nSections; i++){{

            theta = pi * i / nSections;
            x = r*sin(theta);
            y = r*cos(theta);
            sum = sum + (fx * sin(theta) + fy * cos(theta));

        }}

        ff << t << "; " << sum / nSections * 2*pi*r << endl;   

    }}
"""

OPEN_SAVE_STREAM = """
ofstream ff("{save_path}");
"""


class ProblemDefinitionOptions:

    def __init__(

        self, grad_rock, DT, time_range, time_steps,
        vx=0., vy=0., pe=0., save_path=None, n_save=200,
        n_plot=0, n_sections_save=1000,
        log_based_time = True

    ):

        self.grad_rock = grad_rock
        self.DT = DT
        self.vx = vx
        self.vy = vy
        self.pe = pe

        self.max_time = max(time_range)
        self.min_time = min(time_range)
        self.time_steps = time_steps
        self.log_based_time = log_based_time

        self.n_sections_save = n_sections_save
        self.save_path = save_path
        self.n_save = n_save
        self.n_plot = n_plot

    def edp_text(self):

        if self.log_based_time:

            delta_t_dict = LOG_DELTAT

        else:

            delta_t_dict = STD_DELTAT

        return_script = TIME_ITERATION_PARAM.format(

            max_time=self.max_time, min_time=self.min_time, time_steps=self.time_steps,
            delta_t_definition=delta_t_dict["Definition"]

        )

        return_script += PROBLEM_DEF_PARAM.format(

            grad_rock=self.grad_rock, DT=self.DT,
            vx=self.vx, vy=self.vy, Pe=self.pe

        )

        if not (self.vx == 0. and self.vy == 0.):
            return_script += SPEED_ADIMENSIONALIZATION

        else:
            return_script += EMPTY_FIELD_INIT

        return_script += PROBLEM_DEFINITION

        if self.n_plot > 0:
            plot_every_init_txt = PLOT_EVERY_INIT.format(n_plot=self.n_plot)
            display_plt_txt = DISPLAY_PLOTS

        else:
            plot_every_init_txt = """"""
            display_plt_txt = """"""

        if self.save_path is not None:
            open_save_stream = OPEN_SAVE_STREAM.format(save_path=self.save_path.replace("\\", "\\\\"))
            save_output = SAVE_OUTPUT

        else:
            open_save_stream = """"""
            save_output = """"""

        return_script += SOLUTION_AND_SAVING.format(

            plot_every_initialization=plot_every_init_txt,
            n_sections_save=self.n_sections_save, n_save=self.n_save,
            save_path=self.save_path.replace("\\", "\\\\"), display_plot=display_plt_txt,
            save_output=save_output, open_save_stream=open_save_stream,
            delta_t_implementation=delta_t_dict["Implementation"]

        )

        return return_script