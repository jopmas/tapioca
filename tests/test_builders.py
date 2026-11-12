"""Tests of the generic builders (tapioca.inputs: MandyocBuilder, VelocityFieldBuilder, TemperatureFieldBuilder and
GenericTimeConditions)."""

import os

import numpy as np
import pytest

import tapioca as tp
from tapioca.inputs import GenericTimeConditions, InputError, MandyocBuilder, TemperatureFieldBuilder, VelocityFieldBuilder
from tapioca.inputs.svg import AIR, DRY_OLIVINE, WET_OLIVINE, WET_QUARTZ

from conftest import quiet

Lx, Lz, Nx, Nz = 1000e3, 300e3, 101, 31  # 10 km x 10 km


@pytest.fixture
def builder(tmp_path):
    '''
    Air (40 km), crust (35 km), lithospheric mantle (thinner on the right) and asthenosphere.
    Each interface is the top of the layer whose material it carries.
    '''
    b = MandyocBuilder(str(tmp_path / "scenario"), Nx=Nx, Nz=Nz, Lx=Lx, Lz=Lz, verbose=False)
    b.create_interface(-1, 0.0, "air", **AIR.params)
    b.create_interface(0, np.where(b.x < Lx / 2, -160e3, -120e3), "asthenosphere", **WET_OLIVINE.params)
    b.create_interface(1, -75e3, "mantle", **DRY_OLIVINE.params)
    b.create_interface(2, -40e3, "crust", **WET_QUARTZ.params)
    b.create_materials_fields()
    return b


# MandyocBuilder ################################################################

def test_builder_grid(builder):
    assert builder.x[0] == 0 and builder.x[-1] == Lx and len(builder.x) == Nx
    assert builder.z[0] == -Lz and builder.z[-1] == 0 and len(builder.z) == Nz


def test_export_interfaces_writes_in_the_builder_folder(builder):
    builder.export_interfaces(export="none")
    assert os.path.isfile(os.path.join(builder.path, "interfaces.txt"))


def test_interface_needs_a_density(builder):
    with pytest.raises(ValueError, match="density"):
        builder.create_interface(3, -50e3, "no_density", C=1.0)


def test_material_fields_follow_the_interfaces(builder):
    rho = builder.DTree.fields.rho
    left, right = rho.isel(x=0), rho.isel(x=-1)
    at = lambda column, z: float(column.sel(z=z))
    assert at(left, -10e3) == AIR.params["rho"]            # sticky air
    assert at(left, -60e3) == WET_QUARTZ.params["rho"]     # crust
    assert at(left, -150e3) == DRY_OLIVINE.params["rho"]   # lithospheric mantle on the left
    assert at(right, -150e3) == WET_OLIVINE.params["rho"]  # asthenosphere below the thinner mantle on the right
    assert at(left, -250e3) == WET_OLIVINE.params["rho"]   # asthenosphere
    for key in tp.MATERIAL_PARAMETERS:
        assert key in builder.DTree.fields and not (builder.DTree.fields[key] == -999).any()


def test_export_interfaces_with_all_the_parameters(builder):
    builder.export_interfaces(export="none")
    lines = open(f"{builder.path}/interfaces.txt").read().splitlines()
    params = [line.split()[0] for line in lines[:len(tp.MATERIAL_PARAMETERS)]]
    assert params == ["C", "rho", "H", "A", "n", "Q", "V", "k", "weakening_seed", "cohesion_min", "cohesion_max",
                      "friction_angle_min", "friction_angle_max"]
    assert all(len(line.split()) == 1 + 4 for line in lines[:len(params)])  # one value per layer
    data = np.loadtxt(lines[len(params):])
    assert data.shape == (Nx, 3)  # one column per interface
    np.testing.assert_array_equal(data[:, 0], np.where(builder.x < Lx / 2, -160e3, -120e3))  # top of the asthenosphere
    np.testing.assert_array_equal(data[:, 2], -40e3)  # top of the crust


def test_params_keep_the_grid_of_the_builder(builder):
    builder.set_params(params_str="nx = 5\ntime_max = 1.0e6\nsolver = direct")
    assert builder.params["nx"] == Nx and builder.params["n_interfaces"] == 3
    assert builder.params["time_max"] == 1.0e6
    builder.export_params()
    text = open(f"{builder.path}/param.txt").read()
    assert "solver" in text and f"nx" in text and "dimensions" in text


def test_export_fields(builder):
    builder.create_temperature_field(_basic_temperature(builder))
    builder.export_field("temperature")
    T = np.loadtxt(f"{builder.path}/input_temperature_0.txt")
    assert T.shape == (Nx * Nz,)


# VelocityFieldBuilder ############################################################

def test_velocity_builder_conserves_the_flux(builder):
    with quiet():
        v = VelocityFieldBuilder(builder, units="cm/y")
        v.set_region("right", [-120e3, -40e3], vx=1.0)
        assert v.integrate_normal_components() > 0
        v.conservate_veloc("bot")
        assert abs(v.integrate_normal_components()) < 1e-12
        builder.create_velocity_field(v)
    vx = builder.DTree.fields.vx
    assert float(vx.sel(x=Lx, z=-80e3)) == pytest.approx(0.01 / tp.SEC_PER_YEAR)  # 1 cm/yr in m/s
    assert float(abs(vx.isel(x=slice(1, -1), z=slice(1, -1))).max()) == 0  # only the boundaries


@pytest.mark.parametrize("regime, moving_boundary", [("extension", "right"), ("extension", "left"),
                                                     ("extension", "symmetric"), ("compression", "right"),
                                                     ("compression", "left"), ("compression", "symmetric")])
def test_kinematic_regimes_conserve_the_flux(builder, regime, moving_boundary):
    with quiet():
        v = VelocityFieldBuilder(builder, units="cm/y")
        v.set_kinematic_regime(regime, 1.0, moving_boundary=moving_boundary, h_const=140e3)
        assert abs(v.integrate_normal_components()) < 1e-15
    left, right = v.velocs["left"].vx, v.velocs["right"].vx
    one = 0.01 / tp.SEC_PER_YEAR  # 1 cm/yr
    lithosphere = dict(z=-100e3)  # constant velocity (surface at -40 km, h_const = 140 km)
    air = dict(z=0.0)  # more than air_gap (10 km) above the surface: no velocity
    assert float(left.sel(**air)) == 0 and float(right.sel(**air)) == 0
    if moving_boundary == "symmetric":
        np.testing.assert_allclose(right.values, -left.values)
        sign = 1 if regime == "extension" else -1  # outwards or inwards
        assert float(right.sel(**lithosphere)) == pytest.approx(sign * one / 2)
    elif regime == "extension":
        moving = right if moving_boundary == "right" else left
        sign = 1 if moving_boundary == "right" else -1
        assert float(moving.sel(**lithosphere)) == pytest.approx(sign * one)  # outwards
    else:
        moving, fixed = (right, left) if moving_boundary == "right" else (left, right)
        sign = -1 if moving_boundary == "right" else 1
        assert float(moving.sel(**lithosphere)) == pytest.approx(sign * one)  # inwards
        assert np.all(fixed.values == 0)


def test_linear_velocity(builder):
    with quiet():
        v = VelocityFieldBuilder(builder, units="m/s")
        v.linear_velocity("left", "vx", clims=[-Lz, 0], vlims=[1.0, 0.0])
    np.testing.assert_allclose(v.velocs["left"].vx.values, -builder.z / Lz, atol=1e-12)


# TemperatureFieldBuilder ##########################################################

ADIABAT = lambda depth, t_pot=1350.0: t_pot / np.exp(-10.0 * 3.28e-5 * depth / 1250.0)  # depth below the surface


def _basic_temperature(builder, t_pot=1350.0):
    """LAB: top of the asthenosphere (ID 0); surface: top of the crust (ID 2, the default)."""
    t = TemperatureFieldBuilder(builder, LAB_interface=0)
    t.apply_linear_temperature(t_pot=t_pot)
    return t


def test_temperature_builder_lithosphere(builder):
    t = TemperatureFieldBuilder(builder, LAB_interface=0)
    np.testing.assert_array_equal(t.surface_interface.values, -40e3)
    np.testing.assert_array_equal(t.lithosphere_thickness.values, np.where(builder.x < Lx / 2, 120e3, 80e3))


def test_linear_temperature(builder):
    T = _basic_temperature(builder).temperature
    assert float(T.min()) == 0.0
    assert float(T.sel(z=0).max()) == 0.0  # sticky air
    deep = T.sel(z=-250e3)  # asthenosphere: adiabat
    np.testing.assert_allclose(deep.values, ADIABAT(210e3))


def test_steady_columns_follow_the_lithosphere(builder):
    t = _basic_temperature(builder).solve_heat_1D()
    T = t.temperature
    np.testing.assert_allclose(T.sel(z=-40e3).values, 0.0, atol=1e-9)  # surface
    assert float(T.sel(z=-100e3).isel(x=-1)) > float(T.sel(z=-100e3).isel(x=0))  # thinner lithosphere is hotter
    builder.create_temperature_field(t)
    assert builder.DTree.fields.temperature.attrs["t_potential"] == 1350.0


def test_steady_2D_conducts_sideways(builder):
    T1 = _basic_temperature(builder).solve_heat_1D().temperature
    T2 = _basic_temperature(builder)
    T2.solve_heat_2D()
    T2 = T2.temperature
    np.testing.assert_allclose(T2.sel(z=-250e3).values, ADIABAT(210e3))  # asthenosphere kept
    np.testing.assert_allclose(T2.sel(z=0).values, 0.0)  # sticky air kept
    # far from the step of the LAB both solutions agree; at the step the cold side is heated
    far = T2.sel(z=-100e3).isel(x=0)
    assert float(far) == pytest.approx(float(T1.sel(z=-100e3).isel(x=0)), abs=1.0)
    near = dict(z=-100e3, x=Lx / 2 - 2e4)
    assert float(T2.sel(**near)) > float(T1.sel(**near)) + 10.0


def test_solve_heat_column_propagates_one_column(builder):
    t = _basic_temperature(builder)
    t.solve_heat_column(x_inx=0)
    np.testing.assert_array_equal(t.temperature.isel(x=0).values, t.temperature.isel(x=-1).values)


def test_yield_strength_envelope(builder):
    builder.create_temperature_field(_basic_temperature(builder))
    builder.calculate_YSE(sr=1e-15)
    fields = builder.DTree.fields
    for name in ["litho_P", "plastic_YSE", "viscous_YSE", "YSE", "YSE_soft"]:
        assert name in fields
    assert (fields.litho_P >= 0).all()
    assert (fields.YSE <= fields.plastic_YSE).all() and (fields.YSE <= fields.viscous_YSE).all()
    assert (fields.YSE_soft <= fields.YSE + 1e-6).all()  # softening lowers the strength


# Time conditions ###################################################################


@pytest.mark.parametrize("param, conditions, kwargs, text", [
    ("scale_bcv", {75.0: 0.01, 25.0: -1.0}, {}, "2\n25.0 -1.0\n75.0 0.01\n"),  # sorted by time
    ("multi_veloc", {10.0: None, 20.0: None}, {}, "2\n10.0\n20.0\n"),
    ("sedimentation_rate", {0: 25, 20: 50.5}, {}, "2\n0.00 25\n20.00 50.5\n"),
    ("sediment_layers", [(0, 5), (1, 6)], {}, "2 5\n0.00 5\n1.00 6\n"),
    ("base_level", {0: -500, 70: -500}, {}, "2\n0.00 -500\n70.00 -500\n"),
    ("base_level", {0: -500, 70: -500}, {"header": None}, "0.00 -500\n70.00 -500\n"),
    ("SP_vc", {0: (5, 25, 35000.0), 1: (6, 25, 35000.0)}, {}, "0.00 5 25.00 35000.0\n1.00 6 25.00 35000.0\n"),
    ("my_condition", {1: (2, 3)}, {"formats": ("{:g}", "{:d}", "{:.1f}")}, "1\n1 2 3.0\n"),
    ("my_condition", {1: 2}, {"header": lambda lines: f"# {len(lines)} line"}, "# 1 line\n1.0 2\n"),
])
def test_generic_time_conditions_text(param, conditions, kwargs, text):
    assert GenericTimeConditions(param, conditions, **kwargs).text() == text


def test_generic_time_conditions_file_names():
    assert GenericTimeConditions("scale_bcv", {0: 1}).fname == "scale_bcv.txt"
    assert GenericTimeConditions("my_condition", {0: 1}).fname == "my_condition.txt"
    assert GenericTimeConditions("base_level", {0: 1}, fname="sea.txt").fname == "sea.txt"


@pytest.mark.parametrize("conditions, kwargs, message", [
    ({}, {}, "no time condition"),
    ({0: (1, 2), 1: 3}, {}, "same number of values"),
    ([1, 2], {}, "time: value"),
    ({0: 1}, {"formats": ("{}",)}, "1 formats for 2 columns"),
    ({0: 1}, {"header": "lines"}, "header must be"),
])
def test_generic_time_conditions_errors(conditions, kwargs, message):
    with pytest.raises(InputError, match=message):
        GenericTimeConditions("my_condition", conditions, **kwargs)


def test_builder_time_conditions(tmp_path):
    builder = MandyocBuilder(str(tmp_path), Nx=3, Nz=3, Lx=10.0, Lz=10.0)
    builder.set_time_conditions("sedimentation_rate", {0: 25})
    builder.set_time_conditions("sedimentation_rate", {0: 30})  # replaces the file set before
    builder.set_time_conditions("base_level", {0: -100})
    written = builder.export_time_conditions()
    assert sorted(os.path.basename(f) for f in written) == ["base_level.txt", "sedimentation_rate.txt"]
    assert (tmp_path / "sedimentation_rate.txt").read_text() == "1\n0.00 30\n"


def test_generic_time_conditions_table():
    """The conditions are stored as a table: the column time and the named columns of the values, sorted by time."""
    import pandas as pd
    sp_vc = GenericTimeConditions("SP_vc", {1: (6, 25, 35000.0), 0: (5, 25, 35000.0)})
    assert list(sp_vc.table.columns) == ["time", "layer", "rate", "Ld"]
    assert sp_vc.table["time"].tolist() == [0.0, 1.0] and sp_vc.table["layer"].tolist() == [5, 6]
    # the same conditions from a DataFrame (columns in any order), the time as index, or a dict of columns
    df = pd.DataFrame({"rate": [25, 25], "layer": [6, 5], "Ld": [35000.0, 35000.0], "time": [1.0, 0.0]})
    assert GenericTimeConditions("SP_vc", df).text() == sp_vc.text()
    assert GenericTimeConditions("SP_vc", df.set_index("time")).text() == sp_vc.text()
    assert GenericTimeConditions("SP_vc", df.to_dict(orient="list")).text() == sp_vc.text()
    # a DataFrame with other names takes the names of the file
    rate = GenericTimeConditions("sedimentation_rate", pd.DataFrame({"time": [0, 20], "Qs": [25, 50]}))
    assert list(rate.table.columns) == ["time", "rate"] and rate.text() == "2\n0.00 25\n20.00 50\n"
    # other conditions: columns given or value_1, value_2, ...
    assert list(GenericTimeConditions("mine", {0: (1, 2)}).table.columns) == ["time", "value_1", "value_2"]
    assert list(GenericTimeConditions("mine", {0: 1}, columns=["v"]).table.columns) == ["time", "v"]


@pytest.mark.parametrize("param, text", [("scale_bcv", "4\n25.0 -1.0\n75.0 0.01\n85.0 -100.0\n145.0 0.01\n"),
                                         ("sediment_layers", "2 5\n0.00 5\n1.00 6\n"),
                                         ("SP_vc", "0.00 5 25.00 35000.0\n2.50 5 40.00 35000.0\n"),
                                         ("multi_veloc", "2\n10.0\n20.0\n")])
def test_generic_time_conditions_read(tmp_path, param, text):
    """A file of Mandyoc read back gives the same table and the same text."""
    path = tmp_path / "input.txt"
    path.write_text(text)
    conditions = GenericTimeConditions.read(param, str(path))
    assert conditions.text() == text and conditions.fname == "input.txt"


def test_generic_time_conditions_table_errors():
    import pandas as pd
    with pytest.raises(InputError, match="named 'time'"):
        GenericTimeConditions("sedimentation_rate", pd.DataFrame({"t": [0], "rate": [1]}))
    with pytest.raises(InputError, match="columns are"):
        GenericTimeConditions("SP_vc", {0: (5, 25)})
