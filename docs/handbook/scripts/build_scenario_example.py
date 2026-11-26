"""Scenario built with the generic builders of tapioca.inputs (example of the handbook).

    tapioca-svg examples/svg_scenarios/rifting/rifting.svg -o template --no-scripts --no-zip   # param.txt with every key
    python build_scenario_example.py                                                             # writes my_scenario/
"""
import numpy as np
from tapioca.inputs import MandyocBuilder, TemperatureFieldBuilder, VelocityFieldBuilder
from tapioca.inputs.svg import AIR, WET_OLIVINE, DRY_OLIVINE, WET_QUARTZ
from tapioca._aux_functions import read_params

Lx, Lz, Nx, Nz = 1000e3, 300e3, 201, 61                  # 5 km x 5 km
b = MandyocBuilder("my_scenario", Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, verbose=False)
# interface i = top of layer i (z positive upwards, 0 at the top); -1 is the sticky air
b.create_interface(-1, 0.0, "air", **AIR.params)
b.create_interface(0, np.where(b.x < Lx / 2, -200e3, -160e3), "asthenosphere", **WET_OLIVINE.params)
b.create_interface(1, -75e3, "mantle", **DRY_OLIVINE.params)
b.create_interface(2, -40e3, "crust", **{**WET_QUARTZ.params, "H": 9e-10})
b.create_materials_fields()

T = TemperatureFieldBuilder(b, LAB_interface=0, surface_interface=2)   # LAB: top of the asthenosphere
T.apply_linear_temperature(t_pot=1350.0)
T.solve_heat_2D()                                                    # steady-state conduction in the lithosphere
b.create_temperature_field(T)

v = VelocityFieldBuilder(b, units="cm/y")
v.set_kinematic_regime("extension", 1.0, moving_boundary="right", h_const=220e3, surface=2)
b.create_velocity_field(v)

b.export_interfaces(export="none")
params = read_params("template/param.txt")              # a param.txt of tapioca-svg, with every key of Mandyoc
params.update(step_max=5000, time_max=20.0e6, step_print=50)
b.set_params(**params)                                  # nx, nz, lx, lz and n_interfaces come from the builder
b.export_params()
b.export_field("temperature", header="T1\nT2\nT3\nT4")
b.export_field("velocity", header="v1\nv2\nv3\nv4")
