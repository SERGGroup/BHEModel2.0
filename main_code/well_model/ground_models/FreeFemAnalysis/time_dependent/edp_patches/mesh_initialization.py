import os.path

GEOMETRY_PARAM = """
// Geometry Parameters
real ratioL={ratio_L}, ratioH={ratio_H};
real r=1., L=r*ratioL, H=r*ratioH;
"""

MESH_PARAM = """
// Base Mesh Parameters
int pointsH={n_points}, pointsCircle={n_points_circle};
real refIndex={ref_index};
int pointsL = pointsH*L/H;

"""

VISUALIZATION_PARAM = """
// Visualization Mesh Parameters
real ratioR={graph_r_ratio};

"""

INDEX_DEFINITION = """
int Cside=99, Ctop=98, Cbottom=97,CMeshExt=96, CMeshInt=95, Ctube=94;

// Define mesh boundary
border C0(t=0, 2*pi){{x=r*sin(t); y=r*cos(t); label=Ctube;}}
border C10(t=0, 2*pi){x=ratioR*r*sin(t); y=ratioR*r*cos(t); label=CMeshExt;}
border C1(t=L, -L){{x=t; y=H; label=Ctop;}}
border C2(t=H, -H){{x=-L; y=t; label=Cside;}}
border C3(t=-L, L){{x=t; y=-H; label=Cbottom;}}
border C4(t=-H, H){{x=L; y=t; label=Cside;}}
"""

MESH_INITIALIZATION = """
mesh Th = buildmesh(C0(pointsCircle)+C10(pointsCircle)+C1(pointsL)+C2(pointsH)+C3(pointsL)+C4(pointsH));
{refine_mesh}"""

SAVE_MESH = """
savemesh(Th,"{mesh_path}");

"""

RETRIEVE_MESH = """
// The triangulated domain Th is on the left side of its boundary
mesh Th = readmesh("{mesh_path}");

"""

VISUALIZATION_MESH = """
// Mesh for visualization purposes
border C01(t=2*pi, 0){{x=ratioR*r*sin(t); y=ratioR*r*cos(t); label=CMeshExt;}}
mesh ThPlot = buildmesh(C0(10)+C01(100));
"""

MESH_REFINEMENT = """
Th = splitmesh(Th, 2 + int(refIndex * r / sqrt(x*x + y*y)));
"""

PLOT_MESH = """
plot(Th, wait=1);
"""

class MeshOptions:

    def __init__(

        self, n_points, n_points_circle, graph_r_ratio,
        ratio_L=1500., ratio_H=500., ref_index=1.,
        mesh_path=None, retrieve_mesh=False,
        add_visualization_mesh=False,
        show_mesh=False

    ):

        self.n_points = n_points
        self.n_points_circle = n_points_circle
        self.graph_r_ratio = graph_r_ratio
        self.retrieve_mesh = retrieve_mesh
        self.ref_index = ref_index

        self.add_visualization_mesh = add_visualization_mesh
        self.mesh_path = mesh_path
        self.show_mesh = show_mesh
        self.ratio_L = ratio_L
        self.ratio_H = ratio_H

        if self.mesh_path is not None:
            self.mesh_path_name = os.path.basename(self.mesh_path)

        else:
            self.mesh_path_name = "mesh_out.mesh"

    def edp_text(self):

        return_script = GEOMETRY_PARAM.format(ratio_L=self.ratio_L, ratio_H=self.ratio_H)
        return_script += MESH_PARAM.format(

            n_points=self.n_points, n_points_circle=self.n_points_circle,
            ref_index=self.ref_index

        )
        return_script += VISUALIZATION_PARAM.format(graph_r_ratio=self.graph_r_ratio)
        return_script += INDEX_DEFINITION

        if self.retrieve_mesh and os.path.isfile(self.mesh_path):
            return_script += RETRIEVE_MESH.format(mesh_path=self.mesh_path.replace("\\", "\\\\"))

        else:

            if self.ref_index > 1:
                refine_mesh_txt = MESH_REFINEMENT

            else:
                refine_mesh_txt = """"""

            return_script += MESH_INITIALIZATION.format(refine_mesh=refine_mesh_txt)

            if self.mesh_path is not None:
                return_script += SAVE_MESH.format(mesh_path=self.mesh_path.replace("\\", "\\\\"))

        if self.show_mesh:
            return_script += PLOT_MESH

        if self.add_visualization_mesh:
            return_script += VISUALIZATION_MESH

        return return_script
