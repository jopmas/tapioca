"""Tests of the SVG scenarios (tapioca.inputs.svg.SVGScenario) and of the command line tools."""

import re
import shutil
import sys
import zipfile

import numpy as np
import pytest

import tapioca as tp
from tapioca.inputs import MandyocBuilder, TemperatureFieldBuilder, VelocityFieldBuilder
from tapioca.inputs.svg import SVGInputError, SVGScenario
from tapioca.inputs.svg import cli, reader

from tapioca import read_preset_temperature

from conftest import COARSE, CONTINENTAL_RIFT, SVG_SCENARIOS, build_scenario, quiet, svg_path


@pytest.fixture(scope="module", params=list(SVG_SCENARIOS))
def written(request, tmp_path_factory):
    '''Each example on the coarse grid, written with write_all: (name, scenario, output folder).'''
    name = request.param
    scenario = build_scenario(name)
    folder = tmp_path_factory.mktemp(name)
    with quiet():
        scenario.write_all(str(folder), plots=False, system="linux", ncores=4)
    return name, scenario, folder


def _read_params(folder):
    return tp.evaluate_params(params_str=(folder / "param.txt").read_text())


# Input files ######################################################################

def test_input_files(written):
    name, scenario, folder = written
    expected = ["interfaces.txt", "param.txt", "input_temperature_0.txt", f"{name}.svg",
                "run_linux.sh", "0_organize_outputs.sh", "zipper.sh", f"infos_{folder.name}.txt"]
    if scenario.velocity_from_ascii:
        expected += ["input_velocity_0.txt", "vel_bc.txt", "velz_bc.txt"]
    if scenario.variable_bcv:
        expected.append("scale_bcv.txt")
    for fname in expected:
        assert (folder / fname).is_file(), fname
    assert (folder / "input_velocity_0.txt").exists() == (SVG_SCENARIOS[name][0] != "quiescence")


def test_inputs_zip(written):
    '''write_all zips the input files in inputs_<folder>.zip.'''
    name, scenario, folder = written
    with zipfile.ZipFile(folder / f"inputs_{folder.name}.zip") as zf:
        zipped = sorted(zf.namelist())
        assert zf.testzip() is None
    inputs = sorted(f.name for f in folder.iterdir() if f.is_file() and f.suffix != ".zip")
    assert zipped == inputs


def test_inputs_zip_apart_from_the_zip_of_the_run(tmp_path):
    '''The zip of the inputs does not touch <folder>.zip (inputs and outputs, made by zipper.sh after the run).'''
    scenario = build_scenario("rifting")
    run_zip = tmp_path / f"{tmp_path.name}.zip"
    with zipfile.ZipFile(run_zip, "w") as zf:
        zf.writestr("time/time_0.txt", "output")
    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    with zipfile.ZipFile(run_zip) as zf:
        assert zf.namelist() == ["time/time_0.txt"]
    with zipfile.ZipFile(tmp_path / f"inputs_{tmp_path.name}.zip") as zf:
        assert "param.txt" in zf.namelist() and "time/time_0.txt" not in zf.namelist()
    with quiet():
        scenario.write_all(str(tmp_path / "no_zip"), plots=False, scripts=False, zip_inputs=False)
    assert not list((tmp_path / "no_zip").glob("*.zip"))


def test_run_script_post_processing(written):
    name, scenario, folder = written
    run = (folder / "run_linux.sh").read_text()
    assert re.search(r"@[A-Z_]+@", run) is None
    for command in ["bash 0_organize_outputs.sh", "meshes2netcdf.jl", "particles2netcdf.jl", "plot_snapshots.py",
                    "frames_generator.py", "bash zipper.sh", "unzip -tq", "comm -23"]:
        assert command in run, command
    assert run.index("mpirun") < run.index("meshes2netcdf.jl") < run.index("bash zipper.sh") < run.index("rm -rf")


def test_param_txt_matches_the_scenario(written):
    _, scenario, folder = written
    params = _read_params(folder)
    assert (params["nx"], params["nz"]) == (scenario.Nx, scenario.Nz)
    assert params["n_interfaces"] == len(scenario.layers) - 1
    assert params["velocity_from_ascii"] == str(scenario.velocity_from_ascii)
    assert params["magmatism_Ts0"] == scenario.magmatism_Ts0


def test_interfaces_txt(written):
    _, scenario, folder = written
    lines = (folder / "interfaces.txt").read_text().splitlines()
    n_params = len(tp.MATERIAL_PARAMETERS)
    assert sorted(line.split()[0] for line in lines[:n_params]) == sorted(tp.MATERIAL_PARAMETERS)
    values = {line.split()[0]: [float(v) for v in line.split()[1:]] for line in lines[:n_params]}
    assert values["rho"] == [layer.params["rho"] for layer in scenario.layers]
    assert values["C"] == [layer.params["C"] for layer in scenario.layers]
    seed = [i for i, layer in enumerate(scenario.layers) if layer.layer_label == scenario.seed_label]
    assert all(values["weakening_seed"][i] == (1.0 if i in seed else -1.0) for i in range(len(scenario.layers)))
    depths = np.loadtxt(lines[n_params:])
    assert depths.shape == (scenario.Nx, len(scenario.layers) - 1)
    expected = -np.array([layer.interface for layer in scenario.layers[1:]]).T
    np.testing.assert_allclose(depths, expected, rtol=0, atol=0.05)  # written with one decimal (%.1f)


def test_initial_temperature(written):
    _, scenario, folder = written
    T = np.loadtxt(folder / "input_temperature_0.txt")
    assert T.shape == (scenario.Nx * scenario.Nz,)
    T = T.reshape(scenario.Nz, scenario.Nx)
    np.testing.assert_array_equal(T, scenario.T)
    assert T.min() == 0.0 and T.max() < scenario.TP + 300
    air = scenario.Z < scenario.thickness_sa[None, :]
    assert np.all(T[air] == 0.0)


def test_velocity_field(written):
    name, scenario, folder = written
    regime, moving_boundary = SVG_SCENARIOS[name]
    if regime == "quiescence":
        assert scenario.VX is None and scenario.velocity_builder is None
        return
    with quiet():
        flux = scenario.velocity_builder.integrate_normal_components()
    assert abs(flux) < 1e-15, "the inflow must equal the outflow"
    v = np.loadtxt(folder / "input_velocity_0.txt").reshape(scenario.Nz, scenario.Nx, 2)
    np.testing.assert_array_equal(v[..., 0], scenario.VX)
    np.testing.assert_array_equal(v[..., 1], scenario.VZ)
    left, right = scenario.VX[:, 0], scenario.VX[:, -1]
    interior = scenario.VX[:, 1:-1]
    assert np.all(interior == 0), "velocity only at the boundaries"
    lithosphere = scenario.Z[:, 0] < scenario.thickness_sa[0] + scenario.h_v_const
    v_ms = scenario.velocity * 0.01 / tp.SEC_PER_YEAR  # cm/yr -> m/s
    if moving_boundary == "symmetric":
        np.testing.assert_allclose(right, -left)
        assert right[lithosphere & (scenario.Z[:, 0] > scenario.thickness_sa[0])].max() == pytest.approx(v_ms / 2)
    elif regime == "compression":
        fixed = right if moving_boundary == "left" else left
        assert np.all(fixed == 0), "the opposite boundary is fixed"
    else:
        moving = right if moving_boundary == "right" else left
        assert abs(moving).max() > 0


@pytest.mark.parametrize("quiescence_time, n_changes", [(0, 3), (300, 4)])
def test_scale_bcv_of_accordion(tmp_path, quiescence_time, n_changes):
    scenario = build_scenario("accordion", quiescence_time=quiescence_time)
    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    values = (tmp_path / "scale_bcv.txt").read_text().split()
    assert int(values[0]) == n_changes and len(values) == 1 + 2 * n_changes


# Settings ##########################################################################

def test_builders_of_the_scenario():
    scenario = build_scenario("rifting")
    assert isinstance(scenario.builder, MandyocBuilder)
    with quiet():
        scenario.build_temperature()
        scenario.build_velocity()
    assert isinstance(scenario.temperature_builder, TemperatureFieldBuilder)
    assert isinstance(scenario.velocity_builder, VelocityFieldBuilder)
    assert scenario.yield_strength_envelope(10).shape == (scenario.Nz,)
    for name in ["rho", "temperature", "vx", "vz", "YSE"]:
        assert name in scenario.builder.DTree.fields


def test_velocity_changed_by_script_is_written(tmp_path):
    scenario = build_scenario("rifting")
    with quiet():
        scenario.build_velocity()
        scenario.velocity_builder.set_region("top", [0, scenario.Lx], vz=1e-12, mode="add")
        scenario.builder.create_velocity_field(scenario.velocity_builder)
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    velz = np.loadtxt(tmp_path / "velz_bc.txt")
    assert velz == pytest.approx(scenario.builder.DTree.fields.vz.isel(z=-1).values)


@pytest.mark.parametrize("moving_boundary", ["left", "right", "symmetric"])
def test_moving_boundary_argument_overrides_the_svg(moving_boundary):
    scenario = build_scenario("rifting", moving_boundary=moving_boundary)
    assert scenario.moving_boundary == moving_boundary and scenario.source["moving_boundary"] == "argument"


@pytest.mark.parametrize("moving_boundary", ["left", "right", "symmetric"])
def test_compression_moving_boundaries(moving_boundary):
    scenario = build_scenario("compression", moving_boundary=moving_boundary)
    with quiet():
        scenario.build_velocity()
        assert abs(scenario.velocity_builder.integrate_normal_components()) < 1e-15
    left, right = scenario.VX[:, 0], scenario.VX[:, -1]
    rock = scenario.Z[:, 0] > scenario.thickness_sa[0] + 5.0e3  # below the surface (Z is the depth): max. inwards velocity in the lithosphere
    v_ms = scenario.velocity * 0.01 / tp.SEC_PER_YEAR  # cm/yr -> m/s
    if moving_boundary == "symmetric":
        np.testing.assert_allclose(right, -left)
        assert left[rock].max() == pytest.approx(v_ms / 2)
        litho = scenario.thickness_litho if(scenario.h_v_reference == 'lithosphere') else scenario.thickness_mantle
        assert scenario.h_v_const == pytest.approx(max(litho[0], litho[-1]) + 20.0e3)
    else:
        moving, fixed, inwards = (left, right, 1) if moving_boundary == "left" else (right, left, -1)
        assert np.all(fixed == 0), "the opposite boundary is fixed"
        assert (inwards * moving[rock]).max() == pytest.approx(v_ms)


def test_param_txt_keyword_arguments(tmp_path):
    scenario = build_scenario("quiescence")
    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False, particles_per_element=100)
    assert _read_params(tmp_path)["particles_per_element"] == 100
    with pytest.raises(SVGInputError, match="Unknown parameter"):
        scenario.parameters(not_a_parameter=1)


@pytest.mark.parametrize("n_params", [7, 13])
def test_preset_initial_temperature(tmp_path, n_params):
    """Preliminary scenario with temperature.nc (read with MandyocScen) and an interfaces.txt of 7 or 13 parameters."""
    for fname in ["param.txt", "temperature.nc"]:
        shutil.copy(CONTINENTAL_RIFT / fname, tmp_path / fname)
    lines = (CONTINENTAL_RIFT / "interfaces.txt").read_text().splitlines()
    if n_params == 13:  # the 6 parameters added to interfaces.txt with the 13 of DEFAULT_MATERIAL
        n_layers = len(lines[0].split()) - 1
        lines = lines[:7] + [f"{key} " + " ".join([str(tp.DEFAULT_MATERIAL[key])] * n_layers)
                             for key in tp.MATERIAL_PARAMETERS[7:]] + lines[7:]
    (tmp_path / "interfaces.txt").write_text("\n".join(lines) + "\n")
    depth, profile, sa, litho, source = read_preset_temperature(str(tmp_path))
    assert np.all(np.diff(depth) > 0) and profile.max() > 1300
    assert sa == pytest.approx(40e3) and litho > 0
    assert "interfaces from" in source
    scenario = build_scenario("rifting", preset_initial_temperature=str(tmp_path))
    with quiet():
        T = scenario.build_temperature()
    assert scenario.thickness_litho_preset == litho
    assert T.min() == 0.0 and T.max() == pytest.approx(profile.max(), rel=0.05)


WATER = """<rect id="water" inkscape:label="water" x="0" y="0" width="1600" height="45" style="fill:#4f9fd8;stroke:none">
      <desc id="desc_water">air
C = 1.0
rho = 1000.0
H = 0.0
h = 2 km</desc>
    </rect>
    <rect id="air\""""


@pytest.fixture
def water_svg(tmp_path):
    """The rifting example with a water layer (2 km) between the sticky air and the upper crust."""
    svg = tmp_path / "rifting_water.svg"
    svg.write_text(svg_path("rifting").read_text().replace('<rect id="air"', WATER, 1))
    return svg


def _labels(scenario):
    return [layer.layer_label for layer in scenario.layers]


def test_magmatism_extraction_adds_a_basalt_layer():
    scenario = build_scenario("rifting", magmatism="on", magmatism_extraction="on")
    assert _labels(scenario)[-2:] == ["basalt", "air"]
    assert scenario.magmatic_layer == len(scenario.layers) - 2


def test_basalt_goes_below_the_water(water_svg):
    with quiet():
        scenario = SVGScenario(str(water_svg), dx=10, dz=10, verbose=False, water_layer="water",
                               magmatism="on", magmatism_extraction="on")
    assert _labels(scenario)[-3:] == ["basalt", "water", "air"]
    assert scenario.magmatic_layer == len(scenario.layers) - 3 == scenario.surface_id
    water = scenario.layers[scenario.water_layer]
    np.testing.assert_allclose(scenario.thickness_sa, water.interface)  # the surface is the base of the water
    np.testing.assert_allclose(water.interface - scenario.layers[-1].interface, 2e3)


def test_water_without_basalt(water_svg):
    with quiet():
        scenario = SVGScenario(str(water_svg), dx=10, dz=10, verbose=False, water_layer="water")
    assert _labels(scenario)[-2:] == ["water", "air"] and "basalt" not in _labels(scenario)
    assert scenario.magmatic_layer == scenario.surface_id == len(scenario.layers) - 3


def test_magmatic_layer_given_by_the_user():
    """A layer of the SVG receives the melt: no basalt layer is added."""
    scenario = build_scenario("rifting", magmatism="on", magmatism_extraction="on", magmatic_layer="upper_crust")
    assert "basalt" not in _labels(scenario)
    assert scenario.layers[scenario.magmatic_layer].layer_label == "upper_crust"
    assert build_scenario("rifting", magmatic_layer=2).magmatic_layer == 2


@pytest.mark.parametrize("settings, message", [({"water_layer": "lower_crust"}, "right below the sticky air"),
                                               ({"water_layer": "ocean"}, "no unit labelled 'ocean'"),
                                               ({"magmatic_layer": 99}, "IDs 0 to")])
def test_water_and_magmatic_layer_errors(settings, message):
    with pytest.raises(SVGInputError, match=message):
        build_scenario("rifting", **settings)


def test_missing_svg_is_an_error(tmp_path):
    with pytest.raises(SVGInputError, match="not found"):
        SVGScenario(str(tmp_path / "missing.svg"))


def test_figures(tmp_path):
    scenario = build_scenario("rifting_symmetric")
    with quiet():
        scenario.write_all(str(tmp_path), plots=True, scripts=False)
    for fname in ["numerical_setup.png", "initial_temperature.png"]:
        assert (tmp_path / fname).stat().st_size > 10_000


# LAB ###############################################################################

BODY = """    <rect id="{label}_body" inkscape:label="{label}" x="400" y="{y}" width="200" height="10" style="fill:#000000;stroke:none">
      <desc id="desc_{label}_body">{desc}</desc>
    </rect>
    <rect id="weak_seed\""""


def _rifting_with_body(tmp_path, label, y):
    """The rifting example with a 10 km body of the unit label (with its description) at depth y [km], 400-600 km."""
    text = svg_path("rifting").read_text()
    desc = re.search(rf'inkscape:label="{label}".*?<desc[^>]*>(.*?)</desc>', text, re.S).group(1)
    svg = tmp_path / f"rifting_{label}_body.svg"
    svg.write_text(text.replace('    <rect id="weak_seed"', BODY.format(label=label, y=y, desc=desc), 1))
    with quiet():
        return SVGScenario(str(svg), verbose=False, **COARSE)


## NOTE (Claude): an asthenosphere lens in the lithospheric mantle gives the asthenosphere a second slot, with zero
## thickness out of the lens; the LAB was the top of that slot in every column (60 km of lithosphere everywhere)
def test_lab_with_an_asthenosphere_lens_in_the_mantle(tmp_path):
    scenario = _rifting_with_body(tmp_path, "asthenosphere", 100)
    x = scenario.x / 1e3
    assert len(reader.layer_indices(scenario.layers, "asthenosphere")) == 2
    np.testing.assert_allclose(scenario.thickness_litho[(x > 410) & (x < 590)], 60e3)  #top of the lens
    np.testing.assert_allclose(scenario.thickness_litho[(x < 350) | (x > 650)], 155e3)  #base of the mantle
    with quiet():
        scenario.build_temperature()
    np.testing.assert_allclose(scenario.temperature_builder.lithosphere_thickness.values, scenario.thickness_litho)


def test_lab_above_a_lithospheric_body_in_the_asthenosphere(tmp_path):
    scenario = _rifting_with_body(tmp_path, "lithospheric_mantle", 250)
    assert len(reader.layer_indices(scenario.layers, "asthenosphere")) == 2
    np.testing.assert_allclose(scenario.thickness_litho, 155e3)


## NOTE (Claude): matplotlib fills the holes of a path with several subpaths (a unit around a lens)
@pytest.mark.parametrize("hole_clockwise, fill_rule, hole_inside", [(True, "nonzero", False), (False, "nonzero", True),
                                                                    (False, "evenodd", False)])
def test_path_contains_keeps_the_holes(hole_clockwise, fill_rule, hole_inside):
    from matplotlib.path import Path as MplPath
    outer = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
    hole = [(4, 4), (6, 4), (6, 6), (4, 6), (4, 4)]
    if(hole_clockwise):
        hole = hole[::-1]
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * 3 + [MplPath.CLOSEPOLY]
    mpath = MplPath(np.array(outer + hole, dtype=float), codes * 2)
    inside = reader.path_contains(mpath, np.array([[5.0, 5.0], [2.0, 2.0], [12.0, 5.0]]), fill_rule)
    assert inside.tolist() == [hole_inside, True, False]


# Command line ######################################################################

def test_cli_tapioca_svg(tmp_path):
    with quiet():
        cli.main([str(svg_path("rifting")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots",
                  "--no-scripts", "--regime", "accordion", "--quiescence-time", "10"])
    assert _read_params(tmp_path)["variable_bcv"] == "True"
    assert (tmp_path / "scale_bcv.txt").is_file()


def test_cli_tapioca_svg_no_post_processing(tmp_path):
    with quiet():
        cli.main([str(svg_path("rifting")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots",
                  "--system", "linux", "--no-post-processing"])
    run = (tmp_path / "run_linux.sh").read_text()
    assert "bash 0_organize_outputs.sh" in run and "# bash zipper.sh" in run
    assert "meshes2netcdf.jl" not in run and "rm -rf" not in run


@pytest.mark.parametrize("system, kwargs, lookup", [("mac", {}, True), ("mac", {"petsc_env": "my-petsc"}, True),
                                                     ("mac", {"petsc_dir": "/petsc"}, False),
                                                     ("mac", {"petsc_env": None}, False), ("linux", {}, False)])
def test_run_script_petsc_env_lookup(tmp_path, system, kwargs, lookup):
    '''run_mac.sh looks for the PETSc of the conda environment petsc_env before the one of petsc_dir.'''
    scenario = build_scenario("rifting")
    scenario.write_scripts(str(tmp_path), system=system, ncores=2, **kwargs)
    run = (tmp_path / f"run_{system}.sh").read_text()
    assert re.search(r"@[A-Z_]+@", run) is None
    assert ("conda env list" in run) == lookup
    if(lookup):
        assert f"PETSC_ENV=${{PETSC_ENV:-{kwargs.get('petsc_env', 'petsc-env')}}}" in run
        assert run.index("conda env list") < run.index("PETSC_DIR=${PETSC_DIR:-") < run.index("mpirun -n")


def test_run_script_hypatia(tmp_path, monkeypatch):
    '''run_hypatia.sh is a SLURM job with the settings of HYPATIA and the paths of the cluster.'''
    from tapioca.inputs.svg._variables import HYPATIA
    monkeypatch.setenv("PETSC_DIR", "/petsc/of/this/computer")  # not used on hypatia
    monkeypatch.delenv("TAPIOCA_MAIL_USER", raising=False)
    scenario = build_scenario("rifting")
    folder = tmp_path / "my_scenario"
    scenario.write_scripts(str(folder), system="hypatia", mail_user="name@usp.br")
    run = (folder / "run_hypatia.sh").read_text()
    assert re.search(r"@[A-Z_]+@", run) is None
    lines = run.splitlines()
    sbatch = [line for line in lines if line.startswith("#SBATCH")]
    assert lines[1:1 + len(sbatch)] == sbatch, "the #SBATCH lines come right after the shebang"
    for option in [f"--ntasks={HYPATIA['ncores']}", f"--qos={HYPATIA['qos']}", f"--time {HYPATIA['time']}",
                   "--job-name svg-my_scenario", "--mail-user=name@usp.br", "--mail-type=BEGIN,END,FAIL"]:
        assert any(option in line for line in sbatch), option
    for module in HYPATIA["modules"]:
        assert f"module load {module}" in run
    assert "PETSC_DIR=${PETSC_DIR:-$HOME/opt/petsc}" in run and f"PETSC_ARCH='{HYPATIA['petsc_arch']}'" in run
    assert f"NUMBER_OF_CORES=${{SLURM_NTASKS:-{HYPATIA['ncores']}}}" in run
    assert "TAPIOCA_DIR=${TAPIOCA_DIR:-$HOME/opt/tapioca}" in run
    assert run.index("mpirun -n") < run.index("bash 0_organize_outputs.sh") < run.index("conda activate") \
        < run.index("meshes2netcdf.jl") < run.index("bash zipper.sh")
    assert "caffeinate" not in run and "conda env list" not in run


def test_run_script_hypatia_options(tmp_path, monkeypatch):
    monkeypatch.setenv("TAPIOCA_MAIL_USER", "env@usp.br")
    with quiet():
        cli.main([str(svg_path("rifting")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots",
                  "--system", "hypatia", "--ncores", "24", "--job-name", "test", "--no-post-processing"])
    run = (tmp_path / "run_hypatia.sh").read_text()
    assert "#SBATCH --ntasks=24" in run and "#SBATCH --job-name test" in run
    assert "#SBATCH --mail-user=env@usp.br" in run
    assert "conda activate" not in run and "meshes2netcdf.jl" not in run
    with zipfile.ZipFile(tmp_path / f"inputs_{tmp_path.name}.zip") as zf:
        assert "run_hypatia.sh" in zf.namelist()
    monkeypatch.delenv("TAPIOCA_MAIL_USER")
    build_scenario("rifting").write_scripts(str(tmp_path / "no_mail"), system="hypatia")
    assert "--mail" not in (tmp_path / "no_mail" / "run_hypatia.sh").read_text()


def _thermal_state(folder):
    '''Thermal parameters of param.txt and the initial temperature written in folder.'''
    params = _read_params(folder)
    temperature = np.loadtxt(folder / "input_temperature_0.txt", skiprows=4)
    return [float(params[name]) for name in THERMAL_NAMES], temperature


THERMAL_NAMES = ["heat_capacity", "thermal_expansion_coefficient", "gravity_acceleration"]
THERMAL_VALUES = dict(zip(THERMAL_NAMES, [700.0, 4.0e-5, 9.8]))


def _assert_thermal(scenario):
    '''The thermal settings of the scenario are those of the temperature builder.'''
    tb = scenario.temperature_builder
    assert (tb.c_cap, tb.alpha, -tb.g) == tuple(THERMAL_VALUES.values())
    assert tuple(getattr(scenario, name) for name in THERMAL_NAMES) == tuple(THERMAL_VALUES.values())


@pytest.mark.parametrize("given_to", ["constructor", "write_all", "attributes"])
def test_thermal_settings_agree_in_param_and_temperature(tmp_path, given_to):
    (tmp_path / "default").mkdir()
    scenario = build_scenario("quiescence")
    with quiet():
        scenario.write_all(str(tmp_path / "default"), plots=False, scripts=False)
    assert _thermal_state(tmp_path / "default")[0] == [1250.0, 3.28e-5, 10.0]

    scenario = build_scenario("quiescence", **(THERMAL_VALUES if given_to == "constructor" else {}))
    if given_to == "attributes":
        for name, value in THERMAL_VALUES.items():
            setattr(scenario, name, value)
    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False, **(THERMAL_VALUES if given_to == "write_all" else {}))
    params, temperature = _thermal_state(tmp_path)
    assert params == list(THERMAL_VALUES.values())
    _assert_thermal(scenario)
    assert not np.allclose(temperature, _thermal_state(tmp_path / "default")[1])


def test_ccapacity_is_another_name_of_heat_capacity(tmp_path):
    scenario = build_scenario("quiescence")
    scenario.ccapacity = 700.0
    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    assert _thermal_state(tmp_path)[0][0] == 700.0
    assert scenario.temperature_builder.c_cap == 700.0


def test_thermal_settings_from_the_svg_and_their_priority(tmp_path):
    text = svg_path("quiescence").read_text()
    assert "regime" in text
    svg = tmp_path / "scenario.svg"
    svg.write_text(re.sub(r"(regime\s*=)", "heat_capacity = 700\n\\1", text, count=1))
    with quiet():
        from_svg = SVGScenario(str(svg), verbose=False, dx=COARSE["dx"], dz=COARSE["dz"])
        argument = SVGScenario(str(svg), verbose=False, dx=COARSE["dx"], dz=COARSE["dz"], heat_capacity=800.0)
    assert (from_svg.heat_capacity, from_svg.source["heat_capacity"]) == (700.0, "SVG")
    assert (argument.heat_capacity, argument.source["heat_capacity"]) == (800.0, "argument")
    assert from_svg.source["gravity_acceleration"] == "default"


@pytest.mark.parametrize("settings", [{"heat_capacity": 0}, {"heat_capacity": "abc"}, {"gravity_acceleration": -10},
                                      {"thermal_expansion_coefficient": -1.0e-5}])
def test_thermal_settings_invalid(settings):
    with pytest.raises(SVGInputError):
        build_scenario("quiescence", **settings)


def test_thermal_setting_after_the_temperature_is_an_error(tmp_path):
    scenario = build_scenario("quiescence")
    with quiet():
        scenario.build_temperature()
        with pytest.raises(SVGInputError, match="heat_capacity"):
            scenario.write_all(str(tmp_path), plots=False, scripts=False, heat_capacity=700.0)
        scenario.gravity_acceleration = 9.8
        with pytest.raises(SVGInputError, match="gravity_acceleration"):
            scenario.write_all(str(tmp_path), plots=False, scripts=False)


def test_yield_strength_envelope_uses_the_gravity():
    results = []
    for gravity in (10.0, 5.0):
        scenario = build_scenario("quiescence", gravity_acceleration=gravity)
        with quiet():
            results.append(scenario.yield_strength_envelope(0))
    assert not np.allclose(results[0], results[1], equal_nan=True)


def test_generic_builder_rejects_different_thermal_parameters():
    scenario = build_scenario("quiescence")
    with quiet():
        scenario.build_temperature()
        scenario.builder.set_params(heat_capacity=1250.0)
        with pytest.raises(tp.inputs.checks.InputError, match="heat_capacity"):
            scenario.builder.set_params(heat_capacity=700.0)

    scenario = build_scenario("quiescence")
    with quiet():
        scenario.builder.set_params(heat_capacity=700.0)
        with pytest.raises(tp.inputs.checks.InputError, match="heat_capacity"):
            scenario.build_temperature()


@pytest.mark.parametrize("args", [["--heat-capacity", "700"], ["--heat-capacity=700"], ["--heat_capacity", "700"],
                                  ["--heat_capacity=700"]])
def test_cli_heat_capacity_spellings(tmp_path, args):
    with quiet():
        cli.main([str(svg_path("quiescence")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots",
                  "--no-scripts", "--no-zip", *args])
    assert _thermal_state(tmp_path)[0][0] == 700.0


def test_cli_thermal_options(tmp_path):
    with quiet():
        cli.main([str(svg_path("quiescence")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots", "--no-scripts",
                  "--no-zip", "--thermal_expansion_coefficient=4e-5", "--gravity-acceleration", "9.8"])
    assert _thermal_state(tmp_path)[0] == [1250.0, 4.0e-5, 9.8]


def test_cli_invalid_thermal_option(tmp_path):
    with pytest.raises(SystemExit, match="heat_capacity"):
        cli.main([str(svg_path("quiescence")), "-o", str(tmp_path), "--heat-capacity", "-1"])


def test_cli_tapioca_svg_error_message(tmp_path):
    with pytest.raises(SystemExit, match="Error"):
        cli.main([str(tmp_path / "missing.svg"), "-o", str(tmp_path)])


def test_cli_tapioca_svg_interfaces(tmp_path, monkeypatch):
    out = tmp_path / "layers.txt"
    monkeypatch.setattr(sys, "argv", ["tapioca-svg-interfaces", str(svg_path("rifting")), "--dx", "10", "-o", str(out)])
    with quiet():
        reader.main()
    lines = out.read_text().splitlines()
    assert len(lines) == len(tp.MATERIAL_PARAMETERS) + 161


###############################################################################
# custom_veloc_field: velocity field built by a script
###############################################################################

BASAL_OUTFLOW = '''
from tapioca._aux_functions import get_interface
from tapioca.inputs.generic_input_classes import VelocityFieldBuilder

def build_velocity(scenario):
    surface = get_interface(scenario.builder.DTree.interfaces.ds, scenario.surface_id)
    z_surf = float(surface.values[0])
    v = VelocityFieldBuilder(scenario.builder, units='cm/y')
    v.set_region('left', [z_surf - 100e3, z_surf], vx=scenario.velocity)
    v.conservate_veloc('bot')
    return v

def infos(scenario):
    return ['pushed from the left, outflow through the base']
'''


def custom_script(tmp_path, text=BASAL_OUTFLOW, name="build_veloc.py"):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_custom_veloc_field_builds_the_velocity(tmp_path):
    script = custom_script(tmp_path)
    scenario = build_scenario("compression", moving_boundary="left", custom_veloc_field=str(script))
    assert scenario.source["custom_veloc_field"] == "argument"
    with quiet():
        scenario.build_velocity()
        flux = scenario.velocity_builder.integrate_normal_components()
    assert abs(flux) < 1e-15
    assert scenario.VX[:, -1].max() == 0 and scenario.VX[:, 0].max() > 0  # right fixed, left pushed
    base = scenario.VZ[np.argmax(scenario.z), :]
    assert base[2:-2].max() < 0 and np.allclose(base[base != 0], base[base != 0][0])  # uniform outflow
    assert np.all(scenario.VZ[np.argmin(scenario.z), :] == 0)  # nothing through the top

    folder = tmp_path / "out"
    with quiet():
        scenario.write_all(str(folder), plots=False, scripts=False)
    assert (folder / "build_veloc.py").read_text() == BASAL_OUTFLOW
    with zipfile.ZipFile(folder / "inputs_out.zip") as zf:
        assert "build_veloc.py" in zf.namelist()
    infos = (folder / "infos_out.txt").read_text()
    assert f"custom velocity field: {script}" in infos and "outflow through the base" in infos
    v = np.loadtxt(folder / "input_velocity_0.txt").reshape(scenario.Nz, scenario.Nx, 2)
    np.testing.assert_array_equal(v[..., 1], scenario.VZ)


def test_custom_veloc_field_plot_shows_the_base(tmp_path):
    script = custom_script(tmp_path)
    scenario = build_scenario("compression", moving_boundary="left", custom_veloc_field=str(script))
    with quiet():
        scenario.plot_setup(str(tmp_path / "setup.png"))
    assert (tmp_path / "setup.png").stat().st_size > 0


def test_custom_veloc_field_relative_to_the_svg(tmp_path):
    svg = tmp_path / "compression.svg"
    text = svg_path("compression").read_text()
    text = re.sub(r"(regime\s*=\s*compression)", r"\1\ncustom_veloc_field = build_veloc.py", text, count=1)
    svg.write_text(text)
    custom_script(tmp_path)
    with quiet():
        scenario = SVGScenario(str(svg), verbose=False, **COARSE)
    assert scenario.source["custom_veloc_field"] == "SVG"
    assert scenario.custom_veloc_field == str(tmp_path / "build_veloc.py")


@pytest.mark.parametrize("text, message", [("x = 1\n", "has no function build_velocity"),
                                           ("def build_velocity(scenario):\n    return None\n", "must return a VelocityFieldBuilder")])
def test_custom_veloc_field_bad_script(tmp_path, text, message):
    script = custom_script(tmp_path, text)
    scenario = build_scenario("compression", custom_veloc_field=str(script))
    with pytest.raises(SVGInputError, match=message), quiet():
        scenario.build_velocity()


def test_custom_veloc_field_checks(tmp_path):
    with pytest.raises(SVGInputError, match="not found"):
        build_scenario("compression", custom_veloc_field=str(tmp_path / "missing.py"))
    other = tmp_path / "veloc.txt"
    other.write_text("")
    with pytest.raises(SVGInputError, match="Python script"):
        build_scenario("compression", custom_veloc_field=str(other))
    with pytest.raises(SVGInputError, match="quiescence"):
        build_scenario("compression", regime="quiescence", custom_veloc_field=str(custom_script(tmp_path)))


# Surface processes and sedimentation ###############################################

# files of the Mandyoc with the sedimentation of Theunissen and its water layer (the example uses sp_format = ggciag)
SP_VC = dict(sp_format="SP_vc", water_layer="water")

def test_sedimentation_layers_and_files(tmp_path):
    """The 70 layers of sediments of zero thickness come from a template drawn outside the domain (as 0.25cm_35km_lc0.5)."""
    scenario = build_scenario("rifting_sedimentation", **SP_VC)
    labels = _labels(scenario)
    assert labels[:5] == ["asthenosphere", "lithospheric_mantle", "lower_crust", "upper_crust", "upper_crust_top"]
    assert labels[5:75] == ["sediments"] * 70 and labels[75:] == ["water", "air"]
    assert scenario.sediment_layers == list(range(5, 75)) and scenario.water_layer == 75
    for i in range(5, 76):  # sediments and water have zero thickness
        np.testing.assert_allclose(scenario.layers[i].interface, scenario.layers[-1].interface)
    assert scenario.surface_id == 74
    np.testing.assert_allclose(scenario.thickness_sa, 40e3)
    assert scenario.layers[5].params["rho"] == 2400.0 and scenario.layers[75].params["rho"] == 1000.0

    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    sp_vc = (tmp_path / "SP_vc.txt").read_text().splitlines()
    assert len(sp_vc) == 70
    assert sp_vc[0] == "0.00 5 25.00 35000.0" and sp_vc[-1] == "69.00 74 25.00 35000.0"
    assert (tmp_path / "base_level.txt").read_text() == "0.00 -500\n70.00 -500\n"
    params = (tmp_path / "param.txt").read_text()
    for line in ["n_interfaces", "sp_surface_processes", "sp_mode", "variable_base_level", "water_interface",
                 "strain_sed", "air_threshold_density"]:
        assert re.search(rf"^{line}\s*=", params, re.M), line
    assert re.search(r"^sp_mode\s*= theunissen_sedimentation", params, re.M)
    assert re.search(r"^water_interface\s*= 75", params, re.M)
    assert re.search(r"^n_interfaces\s*= 76", params, re.M)
    names = zipfile.ZipFile(tmp_path / f"inputs_{tmp_path.name}.zip").namelist()
    assert "SP_vc.txt" in names and "base_level.txt" in names


def test_sedimentation_window():
    """ti_sed, tf_sed and dt_sed set the number of layers and the times of SP_vc.txt."""
    scenario = build_scenario("rifting_sedimentation", ti_sed=5, tf_sed=20, dt_sed=2.5, **SP_VC)
    assert scenario.sediment_times == [5.0, 7.5, 10.0, 12.5, 15.0, 17.5]
    assert len(scenario.sediment_layers) == 6
    assert scenario.sp_vc().splitlines()[0] == "5.00 5 25.00 35000.0"
    # tf_sed defaults to time_max
    assert len(build_scenario("rifting_sedimentation", tf_sed=None, time_max=30).sediment_times) == 30


def test_no_sedimentation_without_surface_processes(tmp_path):
    scenario = build_scenario("rifting_sedimentation", sp_surface_processes=False, sp_mode="none")
    assert "sediments" not in _labels(scenario) and scenario.sp_vc() is None
    params = scenario.parameters()
    assert params["sp_surface_processes"] is False and "water_interface" not in params
    # the plain examples are not changed
    plain = build_scenario("rifting")
    assert plain.sediment_layers == [] and "sp_mode" not in plain.parameters()


@pytest.mark.parametrize("settings, message", [
    ({"sp_mode": "none"}, "needs sp_mode"),
    ({"sediment_layer": "basin"}, "no unit labelled 'basin'"),
    ({"sp_mode": "erosion"}, "Unknown sp_mode"),
    ({"tf_sed": 0.5, "dt_sed": 1}, "No layer of sediments"),
    ({"sp_surface_processes": "maybe"}, "True or False")])
def test_sedimentation_errors(settings, message):
    with pytest.raises(SVGInputError, match=message):
        build_scenario("rifting_sedimentation", **settings)


def test_unit_layer_of_a_template():
    layer = reader.unit_layer(str(svg_path("rifting_sedimentation")), "sediments", Nx=11)
    assert layer.title == "Sediments" and layer.color == "#e8c96b" and layer.params["k"] == 7.0e-7
    assert reader.unit_layer(str(svg_path("rifting_sedimentation")), "missing", Nx=11) is None


def test_cli_sedimentation_options(tmp_path):
    with quiet():
        cli.main([str(svg_path("rifting_sedimentation")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots",
                  "--no-scripts", "--no-zip", "--ti-sed", "10", "--dt-sed", "5", "--base-level", "-200",
                  "--sp-format", "SP_vc"])
    assert len((tmp_path / "SP_vc.txt").read_text().splitlines()) == 12
    assert (tmp_path / "base_level.txt").read_text() == "0.00 -200\n70.00 -200\n"


def test_sedimentation_rate_in_steps():
    """The rate changes in steps; a change between two layers adds a line with the active layer."""
    scenario = build_scenario("rifting_sedimentation", tf_sed=6, dt_sed=2,
                              sedimentation_rate="0 Myr: 25, 2 Myr: 50, 3000 kyr: 10", **SP_VC)
    assert scenario.sedimentation_rates == [(0.0, 25.0), (2.0, 50.0), (3.0, 10.0)]
    assert scenario.sp_vc().splitlines() == ["0.00 5 25.00 35000.0", "2.00 6 50.00 35000.0",
                                             "3.00 6 10.00 35000.0", "4.00 7 10.00 35000.0",
                                             "6.00 7 0.00 35000.0"]  # tf_sed = 6 Myr stops the sedimentation
    # the same with a dict, from Python
    same = build_scenario("rifting_sedimentation", tf_sed=6, dt_sed=2, sedimentation_rate={0: 25, 2: 50, 3: 10}, **SP_VC)
    assert same.sp_vc() == scenario.sp_vc()
    # a constant rate is the same in every line
    assert {line.split()[2] for line in build_scenario("rifting_sedimentation", **SP_VC).sp_vc().splitlines()} == {"25.00"}


@pytest.mark.parametrize("rate, message", [("5 Myr: 25", "after ti_sed"), ("0: 25, 0: 30", "same time"),
                                           ("0: 25, 10", "time: rate"), ("0: -1", "at least 0")])
def test_sedimentation_rate_errors(rate, message):
    with pytest.raises(SVGInputError, match=message):
        build_scenario("rifting_sedimentation", sedimentation_rate=rate)


def test_cli_sedimentation_rate_in_steps(tmp_path):
    with quiet():
        cli.main([str(svg_path("rifting_sedimentation")), "-o", str(tmp_path), "--dx", "10", "--dz", "10", "--no-plots",
                  "--no-scripts", "--no-zip", "--sedimentation-rate", "0 Myr: 25, 35 Myr: 5", "--sp-format", "SP_vc"])
    rates = [line.split()[2] for line in (tmp_path / "SP_vc.txt").read_text().splitlines()]
    assert rates[:35] == ["25.00"] * 35 and rates[35:] == ["5.00"] * 35


def test_magmatic_layer_from_a_template():
    """A template given as magmatic_layer goes below the sediments and the water, with its own material."""
    scenario = build_scenario("rifting_sedimentation", magmatism="on", magmatism_extraction="on",
                              magmatic_layer="extrusive_basalt", tf_sed=3, **SP_VC)
    labels = _labels(scenario)
    assert labels[4:] == ["upper_crust_top", "extrusive_basalt", "sediments", "sediments", "sediments", "water", "air"]
    assert scenario.magmatic_layer == 5 and scenario.sediment_layers == [6, 7, 8] and scenario.water_layer == 9
    assert "basalt" not in labels  # the template replaces the basalt of tapioca
    basalt = scenario.layers[5]
    assert basalt.params["rho"] == 2900.0 and basalt.color == "#3b3b3b"
    np.testing.assert_allclose(basalt.interface, scenario.layers[-1].interface)  # zero thickness
    params = scenario.parameters()
    assert params["magmatic_layer"] == 5 and params["water_interface"] == 9
    assert scenario.sp_vc().splitlines()[0] == "0.00 6 25.00 35000.0"


def test_default_basalt_with_sedimentation():
    """Without magmatic_layer, the basalt of tapioca goes below the sediments and the water."""
    scenario = build_scenario("rifting_sedimentation", magmatism="on", magmatism_extraction="on", tf_sed=2, **SP_VC)
    assert _labels(scenario)[4:] == ["upper_crust_top", "basalt", "sediments", "sediments", "water", "air"]
    assert scenario.magmatic_layer == 5 and scenario.water_layer == 8


def test_sp_format_ggciag(tmp_path):
    """The files of the Mandyoc of ggciag: sediment_layers.txt, sedimentation_rate.txt and base_level.txt with a header."""
    with quiet():
        scenario = build_scenario("rifting_sedimentation", sp_format="ggciag", tf_sed=4, dt_sed=2,
                                  sedimentation_rate="0: 25, 1: 50")
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    assert not (tmp_path / "SP_vc.txt").exists()
    assert (tmp_path / "sediment_layers.txt").read_text() == "2 5\n0.00 5\n2.00 6\n"
    # rate 0 at tf_sed stops the sedimentation before time_max
    assert (tmp_path / "sedimentation_rate.txt").read_text() == "3\n0.00 25\n1.00 50\n4.00 0\n"
    assert (tmp_path / "base_level.txt").read_text() == "2\n0.00 -500\n70.00 -500\n"
    params = (tmp_path / "param.txt").read_text()
    assert "water_interface" not in params and re.search(r"^variable_base_level\s*= True", params, re.M)
    names = zipfile.ZipFile(tmp_path / f"inputs_{tmp_path.name}.zip").namelist()
    assert "sediment_layers.txt" in names and "sedimentation_rate.txt" in names


def test_sp_format_ggciag_layers_with_diffusion():
    """The Mandyoc of ggciag reads sediment_layers.txt with any sp_mode, but sedimentation_rate.txt only when it deposits."""
    scenario = build_scenario("rifting_sedimentation", sp_format="ggciag", sp_mode="diffusion", tf_sed=2)
    assert set(scenario.sp_files()) == {"sediment_layers.txt", "base_level.txt"}
    assert build_scenario("rifting_sedimentation", sp_mode="diffusion", **SP_VC).sediment_layers == []  # SP_vc: no layers


def test_tf_sed_before_time_max_stops_the_sedimentation():
    lines = build_scenario("rifting_sedimentation", tf_sed=3, **SP_VC).sp_vc().splitlines()
    assert lines == ["0.00 5 25.00 35000.0", "1.00 6 25.00 35000.0", "2.00 7 25.00 35000.0", "3.00 7 0.00 35000.0"]


def test_default_sp_format_is_ggciag(tmp_path):
    """Without sp_format, the files are those of the Mandyoc of ggciag, and the example has no water layer."""
    scenario = build_scenario("rifting_sedimentation", sp_format=None)
    assert scenario.sp_format == "ggciag" and scenario.water_layer is None
    assert _labels(scenario)[-2:] == ["sediments", "air"] and len(scenario.layers) == 76
    assert set(scenario.sp_files()) == {"sediment_layers.txt", "sedimentation_rate.txt", "base_level.txt"}
    assert "water_interface" not in scenario.parameters()


def test_time_conditions_set_by_script_are_written_and_zipped(tmp_path):
    """A file of time conditions set in the builder by script is written by write_all, with those of the scenario."""
    scenario = build_scenario("rifting_sedimentation")
    scenario.builder.set_time_conditions("multi_veloc", {10.0: None})
    with quiet():
        scenario.write_all(str(tmp_path), plots=False, scripts=False)
    assert (tmp_path / "multi_veloc.txt").read_text() == "1\n10.0\n"
    assert (tmp_path / "sediment_layers.txt").exists()
    names = zipfile.ZipFile(tmp_path / f"inputs_{tmp_path.name}.zip").namelist()
    assert "multi_veloc.txt" in names and "sedimentation_rate.txt" in names
    assert [c.param for c in scenario.time_conditions()] == ["sediment_layers", "sedimentation_rate", "base_level"]
