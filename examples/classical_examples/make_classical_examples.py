'''
Writes the SVGs of the classical examples of Mandyoc (mandyoc/examples, one create_input*.ipynb each) in the folders
next to this script:

    python make_classical_examples.py

Each SVG has the geometry, the resolution and the materials of its notebook, and the parameters of param.txt of the
notebook that differ from the defaults of tapioca (in the description of the setup shape). The velocity field of the
punch is built by punch/build_veloc.py (custom_veloc_field), which is not written here, and the initial temperature
of two_layers_with_variable_conductivity (700 * depth / Lz, as in its notebook) is written to
two_layers_with_variable_conductivity/initial_temperature (preset_initial_temperature, field mode).
'''

import os

import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))

###############################################################################
# Materials
###############################################################################

# Materials of the continental rift of Mandyoc (interfaces.txt of create_input_and_plot.ipynb)
ASTHENOSPHERE = """wet olivine
C = 1.0
rho = 3378.0
H = 0.0
A = 1.393e-14
n = 3.0
Q = 429.0e3
V = 15.0e-6"""
LITHO_MANTLE = """dry olivine
C = {C}
rho = 3354.0
H = 9.0e-12
h = {h} km
A = 2.4168e-15
n = 3.5
Q = 540.0e3
V = 25.0e-6"""
LOWER_CRUST = """wet quartz
C = 1.0
rho = 2800.0
H = 0.8e-6 / 2800.0
h = 20 km
A = 8.574e-28
n = 4.0
Q = 222.0e3
V = 0.0"""
UPPER_CRUST = """wet quartz
C = 1.0
rho = 2700.0
H = 2.5e-6 / 2700.0
h = 20 km
A = 8.574e-28
n = 4.0
Q = 222.0e3
V = 0.0"""
AIR = """air
C = {C}
rho = {rho}
H = 0.0
h = {h} km
A = {A}
n = {n}
Q = 0.0
V = 0.0"""
# Sediments of the layers deposited by Mandyoc (template below the domain), as in svg_scenarios/rifting_sedimentation
SEDIMENTS = """wet quartz
C = 1.0
rho = 2400.0
H = 5.417e-10
A = 8.574e-28
n = 4.0
Q = 222.0e3
V = 0.0
k = 7.0e-7
cohesion_min = 5.0e6
cohesion_max = 10.0e6
friction_angle_min = 5.0
friction_angle_max = 15.0"""
# Materials with constant viscosity (rheology_model = 0): the flow law is not used, so A, n, Q and V are 0 as in the
# notebooks; the first line is only the base of the material
LINEAR = "wet olivine\nC = {C}\nrho = {rho}\nH = 0.0{h}\nA = 0.0\nn = 0.0\nQ = 0.0\nV = 0.0"

COLORS = dict(sediments="#e8c96b", asthenosphere="#cfe2cd", lithospheric_mantle="#9bc29b", lower_crust="#f0d1bc",
              upper_crust="#e49c7c", air="#ffffff", weak_seed="#000000", plume="#d9534f", lid="#8fa9c9",
              lithosphere="#9bc29b", dense_layer="#9bc29b", upper_layer="#e49c7c")

###############################################################################
# Parameters of param.txt of the notebooks that differ from the defaults of tapioca
###############################################################################

RIFT_PARAMS = dict(denok="1.0e-15", Xi_min="1.0e-7", sp_surface_tracking="False", step_max="5000", step_print="10",
                   rheology_model="9")
PARAMS = {
    "continental_rift": RIFT_PARAMS,
    "continental_rift_with_surface_processes": {**RIFT_PARAMS, "sp_surface_tracking": "True"},
    "variable_velocity_field": {**RIFT_PARAMS, "step_max": "10000", "step_print": "20"},
    "model1": {**RIFT_PARAMS, "solver": "iterative", "denok": "1.0e-5", "nondimensionalization": "True"},
    "model2": {**RIFT_PARAMS, "solver": "iterative", "denok": "1.0e-9", "rtol": "1.0e-6", "Xi_min": "1.0e-6",
               "sub_division_time_step": "2.0", "viscosity_max": "1.0e24", "viscosity_min": "1.0e19",
               "nondimensionalization": "True"},
    "Crameri2012_case2": dict(denok="1.0e-15", particles_per_element="600", particles_per_element_x="6",
                              particles_per_element_z="100", Xi_min="1.0e-14", random_initial_strain="0.0",
                              initial_dynamic_range="False", basal_heat="-1.0", step_max="7000", step_print="10",
                              sub_division_time_step="1.0", initial_print_step="2", viscosity_reference="1.0e21",
                              viscosity_max="1.0e23", viscosity_min="1.0e19", temperature_from_ascii="False",
                              sticky_blanket_air="False", temperature_difference="0.0", non_linear_method="off",
                              adiabatic_component="off", bot_tangential_velocity="fixed",
                              left_tangential_velocity="free", right_tangential_velocity="free", rheology_model="0",
                              T_initial="0"),
    "slab_detachment": dict(denok="1.0e-15", Xi_min="1.0e-7", sp_surface_tracking="False", step_max="5000",
                            dt_max="1.0e6", step_print="10", initial_print_step="1", viscosity_reference="3.0e22",
                            viscosity_min="1.0e21", viscosity_mean_method="harmonic", temperature_difference="0.0",
                            adiabatic_component="off", radiogenic_component="off", plasticity="off",
                            rheology_model="9"),
    "vanKeken1997_case1a": dict(particles_per_element_x="0", particles_per_element_z="0", denok="1.0e-15",
                                particles_per_element="1100", particles_perturb_factor="0.0", rtol="1.0e-5",
                                Xi_min="1.0e-14", random_initial_strain="0.0", initial_dynamic_range="False",
                                basal_heat="-1.0", sp_surface_tracking="False", step_max="4000", dt_max="1000.0e6",
                                step_print="20", sub_division_time_step="1.0", viscosity_reference="1.0e21",
                                viscosity_min="1.0e17", viscosity_mean_method="harmonic", temperature_from_ascii="False",
                                geoq_fac="1.0", temperature_difference="0.0", non_linear_method="off",
                                adiabatic_component="off", radiogenic_component="off", top_tangential_velocity="fixed",
                                bot_tangential_velocity="fixed", left_tangential_velocity="free",
                                right_tangential_velocity="free", left_temperature="free", right_temperature="free",
                                rheology_model="0", T_initial="0"),
    "two_layers_with_variable_conductivity": dict(denok="1.0e-15", Xi_min="1.0e-7", sp_surface_tracking="False",
                                                  step_max="20000", dt_max="5.0e3", step_print="1000",
                                                  viscosity_min="1.0e25", non_linear_method="off",
                                                  adiabatic_component="off", radiogenic_component="off",
                                                  left_temperature="free", right_temperature="free",
                                                  rheology_model="9"),
    "punch": dict(denok="1.0e-11", particles_per_element="81", Xi_min="1.0e-12", random_initial_strain="0.0",
                  basal_heat="-1.0", sp_surface_tracking="False", step_max="0", step_print="1",
                  viscosity_reference="1.0e4", viscosity_max="1.0e4", viscosity_min="1.0e-2", bcv_extern="True",
                  temperature_from_ascii="False", geoq_fac="1.0", temperature_difference="0.0",
                  adiabatic_component="off", radiogenic_component="off", rheology_model="70", T_initial="0"),
}
PARAMS["model3"] = PARAMS["model2"]


def params_text(name):
    '''Lines of the parameters of param.txt of the example name for the description of the setup shape.'''
    return "".join(f"\n{key} = {value}" for key, value in PARAMS.get(name, {}).items())

###############################################################################
# Drawing
###############################################################################

class Drawing:
    '''SVG of a scenario: shapes of the units ("lithologies" layer), templates outside the domain and the setup shape.'''

    def __init__(self, Lx, Lz, scale=1.0, unit="1 km"):
        self.Lx, self.Lz, self.scale, self.unit = Lx, Lz, scale, unit  # scale: user units per km
        self.shapes, self.templates = [], []
        self.n = 0

    def _common(self, ident, label, title, desc):
        self.n += 1
        fill = COLORS.get(label, "#808080")
        return (f'id="{ident}" inkscape:label="{label}"', fill,
                f'      <title id="title{self.n}">{title}</title>\n      <desc id="desc{self.n}">{desc}</desc>\n')

    def rect(self, ident, label, title, desc, x, y, w, h):
        s = self.scale
        attrs, fill, inner = self._common(ident, label, title, desc)
        self.shapes.append(f'    <rect {attrs} x="{x*s:g}" y="{y*s:g}" width="{w*s:g}" height="{h*s:g}" '
                           f'style="fill:{fill};stroke:none">\n{inner}    </rect>')

    def circle(self, ident, label, title, desc, cx, cy, r):
        s = self.scale
        attrs, fill, inner = self._common(ident, label, title, desc)
        self.shapes.append(f'    <circle {attrs} cx="{cx*s:g}" cy="{cy*s:g}" r="{r*s:g}" '
                           f'style="fill:{fill};stroke:none">\n{inner}    </circle>')

    def path(self, ident, label, title, desc, points):
        s = self.scale
        attrs, fill, inner = self._common(ident, label, title, desc)
        d = "M " + " L ".join(f"{x*s:.4f},{y*s:.4f}" for x, y in points) + " Z"
        self.shapes.append(f'    <path {attrs} d="{d}" style="fill:{fill};stroke:none">\n{inner}    </path>')

    def template(self, ident, label, title, desc, x, y, w, h):
        '''Shape outside the domain ("templates" layer): only gives the material of a layer of zero thickness.'''
        self.rect(ident, label, title, desc, x, y, w, h)
        self.templates.append(self.shapes.pop())

    def _templates(self):
        if not self.templates:
            return ""
        return ('  <g inkscape:groupmode="layer" id="layer_templates" inkscape:label="templates">\n'
                '    <!-- outside the domain: materials of the layers of zero thickness added by tapioca -->\n'
                + "\n".join(self.templates) + "\n  </g>\n")

    def write(self, fname, comments, setup):
        s = self.scale
        W, H = self.Lx * s, self.Lz * s
        head = "\n".join(f"<!-- {c} -->" for c in comments)
        svg = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
{head}
<!-- 1 user unit = {self.unit}; page = model domain. Written by make_classical_examples.py. -->
<svg
   width="{W:g}mm"
   height="{H:g}mm"
   viewBox="0 0 {W:g} {H:g}"
   version="1.1"
   id="svg1"
   inkscape:version="1.4"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   xmlns="http://www.w3.org/2000/svg">
  <g inkscape:groupmode="layer" id="layer1" inkscape:label="lithologies">
{chr(10).join(self.shapes)}
  </g>
{self._templates()}  <g inkscape:groupmode="layer" id="layer_setup" inkscape:label="domain">
    <rect id="setup" inkscape:label="setup" x="0" y="0" width="{W:g}" height="{H:g}" style="fill:none;stroke:#ff0000;stroke-width:{max(W, H)/1000:g}">
      <title id="title_setup">Model setup</title>
      <desc id="desc_setup">{setup}</desc>
    </rect>
  </g>
</svg>
'''
        os.makedirs(os.path.dirname(fname), exist_ok=True)
        with open(fname, "w") as f:
            f.write(svg)
        print("written", os.path.relpath(fname, OUT))

###############################################################################
# Examples
###############################################################################

def seed_columns(Nx, Lx, x_seed, n_seed=2):
    '''x [m] of the first and last nodes of the seed of the notebooks of Mandyoc (n_seed nodes around x_seed) and dx.'''
    i0 = int(Nx * x_seed // Lx - n_seed // 2)
    i1 = int(Nx * x_seed // Lx + n_seed // 2)
    dx = Lx / (Nx - 1)
    return i0 * dx, (i1 - 1) * dx, dx


NOTE_RIFT = ("Layers (from Mandyoc): sticky air 40 km, upper crust 20 km, lower crust 20 km, lithosphere 130 km; "
             "weak seed = lithospheric mantle with C = 0.1, 6 km thick, 13 km below the Moho, two nodes wide at x = {x} km.")
NOTE_TEMPERATURE = ("Difference from Mandyoc: tapioca computes its own initial temperature (conduction with the LAB at "
                    "1350 oC; Mandyoc: 1300 oC at the LAB, Tp = 1262 oC and 500 Myr of diffusion). The velocity of the "
                    "rifting regime is the one of the notebook: constant down to 150 km below the surface (lithosphere + "
                    "20 km), the return flow below it and no velocity in the top 30 km of the sticky air.")


def continental_rift(name, fname, Nx, Nz, comments, extra_setup="", Lx=1600.0, x_seed=750.0, sediments=False,
                     experiment="regime = rifting\nvelocity = 1 cm/yr\nmoving_boundary = right\ntime_max = 100 Myr"):
    Lz = 300.0
    sa, uc, lc, litho, seed_depth, seed_h = 40.0, 20.0, 20.0, 130.0, 13.0, 6.0
    d = Drawing(Lx, Lz)
    d.rect("asthenosphere", "asthenosphere", "Asthenosphere", ASTHENOSPHERE, 0, 0, Lx, Lz)
    d.rect("lithospheric_mantle", "lithospheric_mantle", "Lithospheric mantle",
           LITHO_MANTLE.format(C=1.0, h=litho - uc - lc), 0, 0, Lx, sa + litho)
    d.rect("lower_crust", "lower_crust", "Lower crust", LOWER_CRUST, 0, 0, Lx, sa + uc + lc)
    d.rect("upper_crust", "upper_crust", "Upper crust", UPPER_CRUST, 0, 0, Lx, sa + uc)
    d.rect("air", "air", "Sticky air", AIR.format(C=1.0, rho=1.0, h=sa, A=1.0e-18, n=1.0), 0, 0, Lx, sa)
    # seed: lithospheric mantle with C = 0.1, n_seed nodes wide, centred 13 km below the Moho, 6 km thick
    xa, xb, dx = (v / 1e3 for v in seed_columns(Nx, Lx * 1e3, x_seed * 1e3))
    zc = sa + uc + lc + seed_depth
    d.rect("weak_seed", "weak_seed", "Weak seed", LITHO_MANTLE.format(C=0.1, h=seed_h),
           xa - dx / 2, zc - seed_h / 2, xb - xa + dx, seed_h)
    if sediments:
        d.template("sediments", "sediments", "Sediments", SEDIMENTS, 0, Lz + 10, 120, 20)
    setup = f"Lx = {Lx:g} km\nLz = {Lz:g} km\nNx = {Nx}\nNz = {Nz}\n{experiment}\nmagmatism = off{extra_setup}{params_text(name)}"
    d.write(fname, comments, setup)


def main():
    # ---------------------------------------------------------------- continental rifts
    continental_rift("continental_rift", f"{OUT}/continental_rift/continental_rift.svg", 1601, 301,
                     ["continental_rift of Mandyoc (examples/continental_rift/create_input_and_plot.ipynb).",
                      NOTE_RIFT.format(x=750), "Extension at 1 cm/yr at the right boundary, 100 Myr.", NOTE_TEMPERATURE])
    continental_rift("continental_rift_with_surface_processes",
                     f"{OUT}/continental_rift_with_surface_processes/continental_rift_with_surface_processes.svg", 801, 151,
                     ["continental_rift_with_surface_processes of Mandyoc (examples/continental_rift_with_surface_processes).",
                      NOTE_RIFT.format(x=750), "Extension at 1 cm/yr at the right boundary, 100 Myr, with surface processes.",
                      "The notebook uses the old sp_mode = 5 (erosion with the precipitation.txt profile, k_sp = 0.0008 m/yr, "
                      "sp_dt = 1.0e5 yr), which the current Mandyoc does not have: here sp_mode = diffusion with sp_d_c = 1.0 "
                      "and sea_level = -1000 m, as in the notebook.",
                      "The Mandyoc of ggciag reads sediment_layers.txt whenever sp_surface_processes is True, so tapioca adds "
                      "one layer of zero thickness of the sediments (template below the domain, the material of "
                      "svg_scenarios/rifting_sedimentation) every dt_sed = 10 Myr.", NOTE_TEMPERATURE],
                     extra_setup="\nsp_surface_processes = True\nsp_mode = diffusion\nsp_d_c = 1.0\nsea_level = -1000.0\n"
                                 "sediment_layer = sediments\ndt_sed = 10 Myr", sediments=True)
    continental_rift("variable_velocity_field", f"{OUT}/variable_velocity_field/variable_velocity_field.svg", 401, 76,
                     ["variable_velocity_field of Mandyoc (examples/variable_velocity_field/create_input_and_plot.ipynb).",
                      NOTE_RIFT.format(x=750),
                      "Accordion with the times of scale_bcv.txt of the notebook: extension at 1 cm/yr at the right boundary "
                      "until 20 Myr (ti_convergence), convergence until 70 Myr (ti_quiescence), then the boundaries stop "
                      "(quiescence_scale = 0) until time_max = 100 Myr, with no second extension.", NOTE_TEMPERATURE],
                     experiment="regime = accordion\nvelocity = 1 cm/yr\nmoving_boundary = right\nti_convergence = 20 Myr\n"
                                "ti_quiescence = 70 Myr\nquiescence_time = 30 Myr\nquiescence_scale = 0\ntime_max = 100 Myr")
    for model, (Nx, Nz, Lx, x_seed, extra) in {
            "model1": (801, 161, 1600.0, 750.0, "2D, the scenario of continental_rift run non-dimensional (nondimensionalization = True)."),
            "model2": (201, 41, 1600.0, 750.0, "3D in Mandyoc (Ly = 60 km, Ny = 9): this is the xz section, the same along y."),
            "model3": (161, 41, 1000.0, 500.0, "3D in Mandyoc (Ly = 1000 km, Ny = 161): this is the xz section, the same along y.")}.items():
        continental_rift(model, f"{OUT}/continental_rift_non_dimensional/{model}/{model}.svg", Nx, Nz,
                         [f"continental_rift_non_dimensional/{model} of Mandyoc.", extra, NOTE_RIFT.format(x=f"{x_seed:g}"),
                          "Extension at 1 cm/yr at the right boundary, 100 Myr.", NOTE_TEMPERATURE],
                         Lx=Lx, x_seed=x_seed)

    # ---------------------------------------------------------------- Crameri et al. (2012), case 2
    Lx, Lz = 2800.0, 850.0
    d = Drawing(Lx, Lz)
    d.rect("asthenosphere", "asthenosphere", "Mantle", LINEAR.format(C=1.0, rho=3300.0, h=""), 0, 0, Lx, Lz)
    d.rect("lid", "lid", "Lid", LINEAR.format(C=100.0, rho=3300.0, h="\nh = 100 km"), 0, 0, Lx, 250.0)
    d.rect("air", "air", "Sticky air", AIR.format(C=0.01, rho=0.0, h=150, A=0.0, n=0.0), 0, 0, Lx, 150.0)
    d.circle("plume", "plume", "Plume", LINEAR.format(C=0.1, rho=3200.0, h=""), Lx / 2, Lz - 300.0, 50.0)
    d.write(f"{OUT}/Crameri2012_case2/Crameri2012_case2.svg",
            ["Crameri2012_case2 of Mandyoc (examples/Crameri2012_case2): rising plume below a lid, Crameri et al. (2012), case 2.",
             "Mantle (C = 1, rho = 3300), plume of 100 km in diameter with its centre 300 km above the base (C = 0.1, "
             "rho = 3200), lid 100 km (C = 100), sticky air 150 km (C = 0.01, rho = 0). Constant viscosity "
             "(rheology_model = 0, viscosity_reference = 1e21 Pa s), isothermal (temperature_difference = 0), free slip at "
             "the sides and no slip at the base, as in param.txt of the notebook."],
            f"Lx = {Lx:g} km\nLz = {Lz:g} km\nNx = 281\nNz = 169\nregime = quiescence\ntime_max = 21 Myr{params_text('Crameri2012_case2')}")

    # ---------------------------------------------------------------- slab detachment
    Lx, Lz = 1000.0, 660.0
    d = Drawing(Lx, Lz)
    d.rect("asthenosphere", "asthenosphere", "Mantle",
           "wet olivine\nC = 1.0\nrho = 3150.0\nH = 0.0\nA = 1.0e-21\nn = 1.0\nQ = 0.0\nV = 0.0", 0, 0, Lx, Lz)
    lith = "dry olivine\nC = 1.0\nrho = 3300.0\nH = 0.0\nA = 1.96438e-47\nn = 4.0\nQ = 0.0\nV = 0.0"
    d.rect("lithosphere", "lithosphere", "Lithosphere", lith, 0, 0, Lx, 80.0)
    d.rect("slab", "lithosphere", "Lithosphere", lith, 459.5, 0, 80.0, 330.0)
    d.write(f"{OUT}/slab_detachment/slab_detachment.svg",
            ["slab_detachment of Mandyoc (examples/slab_detachment): Schmalholz (2011).",
             "Lithosphere 80 km (rho = 3300, n = 4) with a vertical slab 80 km wide and 250 km long at the centre, over a "
             "linear mantle (rho = 3150). No sticky air (sticky_air = off): tapioca adds one of zero thickness at the top.",
             "Isothermal: temperature_difference = 0, so Mandyoc caps the temperature of input_temperature_0.txt at 0 oC, "
             "as the T = 0 of the notebook."],
            f"Lx = {Lx:g} km\nLz = {Lz:g} km\nNx = 1001\nNz = 601\nsticky_air = off\nregime = quiescence\ntime_max = 25 Myr\n"
            f"gravity_acceleration = 9.81\nmin_thickness = 0{params_text('slab_detachment')}")

    # ---------------------------------------------------------------- van Keken et al. (1997), case 1a
    Lx, Lz = 91.42, 100.0
    d = Drawing(Lx, Lz)
    d.rect("asthenosphere", "asthenosphere", "Light layer", LINEAR.format(C=1.0, rho=-1000.0, h=""), 0, 0, Lx, Lz)
    x = np.linspace(0, Lx, 201)
    depth = 0.8 * Lz - 0.02 * Lz * np.cos(np.pi * x / Lx)
    d.path("dense_layer", "dense_layer", "Dense layer", LINEAR.format(C=1.0, rho=0.0, h=""),
           [(0.0, 0.0), (Lx, 0.0)] + list(zip(x[::-1], depth[::-1])))
    d.write(f"{OUT}/vanKeken1997_case1a/vanKeken1997_case1a.svg",
            ["vanKeken1997_case1a of Mandyoc (examples/vanKeken1997_case1a): Rayleigh-Taylor instability, van Keken et al. (1997), case 1a.",
             "Interface at z = -0.8 Lz + 0.02 Lz cos(pi x / Lx): light layer (rho = -1000) below, dense layer (rho = 0) "
             "above, the same viscosity (C = 1). No sticky air (sticky_air = off). Isothermal, constant viscosity, no slip "
             "at the top and the base, as in param.txt of the notebook (time_max = 1.4e11, written in s in the notebook)."],
            f"Lx = {Lx:g} km\nLz = {Lz:g} km\nNx = 81\nNz = 81\nsticky_air = off\nregime = quiescence\ntime_max = 140000 Myr\n"
            f"min_thickness = 0{params_text('vanKeken1997_case1a')}")

    # ---------------------------------------------------------------- two layers with variable conductivity
    name = "two_layers_with_variable_conductivity"
    Lx, Lz, Nx, Nz = 40.0, 40.0, 11, 11
    flow = "wet olivine\nC = 1.0\nrho = 3300.0\nH = 0.0{h}\nA = 1.393e-14\nn = 3.0\nQ = 429.0e3\nV = 15.0e-6\nk = {k}"
    d = Drawing(Lx, Lz)
    d.rect("asthenosphere", "asthenosphere", "Lower layer", flow.format(h="", k="1.0e-6"), 0, 0, Lx, Lz)
    d.rect("upper_layer", "upper_layer", "Upper layer", flow.format(h="\nh = 20 km", k="0.25e-6"), 0, 0, Lx, 20.0)
    d.write(f"{OUT}/{name}/{name}.svg",
            [f"{name} of Mandyoc (examples/{name}).",
             "Two layers of 20 km with the same material but the thermal diffusivity: k = 1e-6 (below) and 0.25e-6 m2/s "
             "(above). No sticky air (sticky_air = off).",
             "Initial temperature of the notebook (700 oC * depth / Lz) in initial_temperature/input_temperature_0.txt "
             "(preset_initial_temperature, field mode); free temperature at the sides, no deformation "
             "(viscosity_min = viscosity_max = 1e25 Pa s)."],
            f"Lx = {Lx:g} km\nLz = {Lz:g} km\nNx = {Nx}\nNz = {Nz}\nsticky_air = off\nregime = quiescence\ntime_max = 100 Myr\n"
            f"min_thickness = 0\npreset_initial_temperature = initial_temperature/input_temperature_0.txt\n"
            f"preset_temperature_mode = field{params_text(name)}")
    folder = f"{OUT}/{name}/initial_temperature"
    os.makedirs(folder, exist_ok=True)
    z = np.linspace(Lz * 1e3, 0, Nz)  # depth, from the bottom to the top (order of Mandyoc)
    T = np.ones((Nz, Nx)) * (700.0 * z / (Lz * 1e3))[:, None]
    np.savetxt(f"{folder}/input_temperature_0.txt", T.ravel(), header="T1\nT2\nT3\nT4")
    with open(f"{folder}/param.txt", "w") as f:  # grid of the field (read by tapioca.read_preset_temperature_field)
        f.write(f"nx = {Nx}\nnz = {Nz}\nlx = {Lx * 1e3:g}\nlz = {Lz * 1e3:g}\n")
    print("written", os.path.relpath(folder, OUT))

    # ---------------------------------------------------------------- punch
    d = Drawing(1.0e-3, 0.5e-3, scale=1.0e6, unit="1 mm")  # 1 m x 0.5 m
    d.rect("asthenosphere", "asthenosphere", "Material", LINEAR.format(C=1.0, rho=0.0, h=""), 0, 0, 1.0e-3, 0.5e-3)
    d.write(f"{OUT}/punch/punch.svg",
            ["punch of Mandyoc (examples/punch): rigid punch on a perfectly plastic material (Prandtl), non-dimensional (1 x 0.5).",
             "One material (rho = 0, g = 0) and no sticky air (sticky_air = off, as the second interface of the notebook, "
             "at the top). The vertical sampling of the SVG is 1 mm (sampling = 0.001 m), smaller than the domain.",
             "The punch is the velocity of build_veloc.py (custom_veloc_field): vz = -1 at the top for |x - 0.5| <= 0.08, "
             "and its bcv_0.txt (bcv_extern = True) frees the other nodes, as in create_input_files_punch.ipynb. The "
             "regime is only the base of the custom velocity field (quiescence has no velocity)."],
            "Lx = 1 m\nLz = 0.5 m\nNx = 401\nNz = 201\nsticky_air = off\nsampling = 0.001 m\nmin_thickness = 0\n"
            "regime = compression\ncustom_veloc_field = build_veloc.py\ntime_max = 130 Myr\ngravity_acceleration = 0"
            f"{params_text('punch')}")


if __name__ == "__main__":
    main()
