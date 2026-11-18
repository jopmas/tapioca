"""
Tests of the classical examples of Mandyoc (examples/classical_examples) and of the settings they need: the times of
the accordion, scenarios without sticky air, the vertical sampling of the SVG, parameters of param.txt in the setup
shape, isothermal scenarios and the files written by the script of custom_veloc_field.
"""

import re
import zipfile

import numpy as np
import pytest

from tapioca.inputs.svg import SVGScenario
from tapioca.inputs.svg.reader import read_setup
from tapioca.inputs.checks import SVGInputError

from conftest import EXAMPLES, quiet, svg_path

CLASSICAL = EXAMPLES / "classical_examples"


def classical_path(name):
    '''Path of the SVG of a classical example (name may be folder/model).'''
    return CLASSICAL / name / f"{name.split('/')[-1]}.svg"


def build_classical(name, **kwargs):
    with quiet():
        return SVGScenario(str(classical_path(name)), verbose=False, **kwargs)


def write_classical(name, folder, **kwargs):
    scenario = build_classical(name, **kwargs)
    with quiet():
        scenario.write_all(str(folder), plots=False, scripts=False)
    return scenario


def read_params(folder):
    params = {}
    for line in (folder / "param.txt").read_text().splitlines():
        line = line.split("#")[0]
        if "=" in line:
            key, value = (s.strip() for s in line.split("=", 1))
            params[key] = value
    return params


def modified_svg(tmp_path, name, old, new):
    '''Copy of the SVG of a classical example with old replaced by new in the text.'''
    text = classical_path(name).read_text()
    assert old in text
    svg = tmp_path / f"{name.split('/')[-1]}.svg"
    svg.write_text(text.replace(old, new, 1))
    return svg

# Examples ##########################################################################

@pytest.mark.parametrize("name", ["Crameri2012_case2", "slab_detachment", "vanKeken1997_case1a",
                                  "two_layers_with_variable_conductivity", "punch", "variable_velocity_field"])
def test_classical_examples_are_built(tmp_path, name):
    write_classical(name, tmp_path)
    for fname in ["interfaces.txt", "param.txt", "input_temperature_0.txt"]:
        assert (tmp_path / fname).is_file()


@pytest.mark.parametrize("name", ["continental_rift", "continental_rift_with_surface_processes",
                                  "continental_rift_non_dimensional/model1", "continental_rift_non_dimensional/model2",
                                  "continental_rift_non_dimensional/model3"])
def test_classical_rifts_are_read(name):
    setup = read_setup(str(classical_path(name)))
    assert setup["experiment"]["regime"] == "rifting"
    assert setup["params"]["rheology_model"] == "9"

# Accordion #########################################################################

def test_accordion_with_the_times_of_variable_velocity_field(tmp_path):
    scenario = write_classical("variable_velocity_field", tmp_path)
    assert (scenario.ti_convergence, scenario.ti_quiescence, scenario.quiescence_scale) == (20.0, 70.0, 0.0)
    assert scenario.time_max == pytest.approx(100.0e6) and scenario.dt_rifting2 == 0.0
    values = (tmp_path / "scale_bcv.txt").read_text().split()
    assert [float(v) for v in values] == [2, 20.0, -1.0, 70.0, 0.0]  # as scale_bcv.txt of the notebook
    assert read_params(tmp_path)["variable_bcv"] == "True"


def test_accordion_times_given_as_arguments():
    with quiet():
        scenario = SVGScenario(str(svg_path("accordion")), verbose=False, dx=10, dz=10, ti_convergence=10,
                               ti_quiescence=30, quiescence_time=5, quiescence_scale=0.1, time_max=50)
    assert scenario._bcv_changes() == [(10.0, -1.0), (30.0, 0.1), (35.0, -10.0), (50.0, 0.01)]


@pytest.mark.parametrize("settings, message", [
    (dict(ti_convergence=50, ti_quiescence=40), "ti_quiescence"),
    (dict(quiescence_time=10, quiescence_scale=0, time_max=200), "quiescence_scale = 0"),
    (dict(quiescence_time=10, time_max=50), "ends before the end of the quiescence")])
def test_accordion_times_errors(settings, message):
    with quiet(), pytest.raises(SVGInputError, match=message):
        SVGScenario(str(svg_path("accordion")), verbose=False, dx=10, dz=10, **settings)

# Sticky air ########################################################################

def test_scenario_without_sticky_air(tmp_path):
    scenario = write_classical("two_layers_with_variable_conductivity", tmp_path)
    air = scenario.layers[-1]
    assert air.rheology_name == "air" and np.all(air.interface == 0.0)
    assert np.all(scenario.thickness_sa == 0.0)
    assert [layer.layer_label for layer in scenario.layers] == ["asthenosphere", "upper_layer", "air"]
    assert read_params(tmp_path)["sticky_blanket_air"] == "False"
    k = [line for line in (tmp_path / "interfaces.txt").read_text().splitlines() if line.startswith("k ")][0]
    assert k.split()[1:3] == ["1e-06", "2.5e-07"]
    # initial temperature of the notebook (preset field): 700 oC * depth / Lz
    T = np.loadtxt(tmp_path / "input_temperature_0.txt").reshape(scenario.Nz, scenario.Nx)
    np.testing.assert_allclose(T[:, 0], 700.0 * scenario.z / scenario.Lz, atol=1.0e-6)


def test_sticky_air_off_with_air_at_the_top_is_an_error():
    with quiet(), pytest.raises(SVGInputError, match="sticky_air = off"):
        SVGScenario(str(svg_path("rifting")), verbose=False, dx=10, dz=10, sticky_air="off")


def test_without_sticky_air_off_the_top_unit_must_be_air(tmp_path):
    svg = modified_svg(tmp_path, "vanKeken1997_case1a", "sticky_air = off\n", "")
    with quiet(), pytest.raises(SVGInputError, match="sticky_air = off"):
        SVGScenario(str(svg), verbose=False)

# Sampling ##########################################################################

def test_sampling_of_the_punch(tmp_path):
    scenario = write_classical("punch", tmp_path)
    assert scenario.sampling == pytest.approx(1.0e-3) and scenario.Lz == pytest.approx(0.5)
    assert scenario.gravity_acceleration == 0.0


def test_sampling_larger_than_the_domain_is_an_error():
    with quiet(), pytest.raises(SVGInputError, match="sampling"):
        SVGScenario(str(classical_path("punch")), verbose=False, sampling=100.0)


def test_sampling_argument_overrides_the_svg():
    with quiet():
        scenario = SVGScenario(str(svg_path("rifting")), verbose=False, dx=10, dz=10, sampling=500.0)
    assert scenario.sampling == 500.0

# Velocity of the punch #############################################################

def test_punch_velocity_and_bcv(tmp_path):
    scenario = write_classical("punch", tmp_path)
    X, Z = scenario.X, scenario.Z
    under = np.isclose(Z, 0.0) & (np.abs(X - 0.5) <= 0.08 + 1.0e-9)
    v = np.loadtxt(tmp_path / "input_velocity_0.txt").reshape(-1, 2)
    vz = v[:, 1].reshape(X.shape)
    assert np.all(vz[under] == -1.0) and np.all(vz[~under] == 0.0) and np.all(v[:, 0] == 0.0)
    bcv = np.loadtxt(tmp_path / "bcv_0.txt").reshape(-1, 2)
    assert np.all(bcv[under.ravel()] == 0.0) and bcv[:, 1].reshape(X.shape)[0].sum() == 0.0
    params = read_params(tmp_path)
    assert (params["bcv_extern"], params["rheology_model"], params["step_max"]) == ("True", "70", "0")
    with zipfile.ZipFile(tmp_path / f"inputs_{tmp_path.name}.zip") as zf:
        assert {"bcv_0.txt", "build_veloc.py"} <= set(zf.namelist())

# Parameters of param.txt in the setup shape ########################################

def test_param_txt_from_the_svg_and_keyword_arguments(tmp_path):
    svg = tmp_path / "rifting.svg"
    svg.write_text(re.sub(r"(regime\s*=)", "rheology_model = 9\nviscosity_min = 1.0e19\n\\1",
                          svg_path("rifting").read_text(), count=1))
    with quiet():
        scenario = SVGScenario(str(svg), verbose=False, dx=10, dz=10)
    assert scenario.svg_params == {"rheology_model": 9, "viscosity_min": 1.0e19}
    params = scenario.parameters(viscosity_min=1.0e20)
    assert params["rheology_model"] == 9 and params["viscosity_min"] == 1.0e20


@pytest.mark.parametrize("line, message", [("nx = 10", "set by the domain"), ("not_a_param = 1", "unknown key"),
                                           ("step_max = abc", "expected int")])
def test_param_txt_from_the_svg_errors(tmp_path, line, message):
    svg = tmp_path / "rifting.svg"
    svg.write_text(re.sub(r"(regime\s*=)", f"{line}\n\\1", svg_path("rifting").read_text(), count=1))
    with quiet(), pytest.raises(SVGInputError, match=message):
        SVGScenario(str(svg), verbose=False, dx=10, dz=10)


@pytest.mark.parametrize("name", ["Crameri2012_case2", "slab_detachment", "vanKeken1997_case1a"])
def test_isothermal_scenarios_start_at_0_oC(tmp_path, name):
    scenario = write_classical(name, tmp_path)
    assert scenario.isothermal and np.all(scenario.T == 0.0)
