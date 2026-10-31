"""Tests of the scenarios that start from the final state of a Mandyoc scenario (tapioca-svg-from-state and
preset_temperature_mode = field of SVGScenario)."""

import re

import numpy as np
import pytest
import xarray as xr

import tapioca as tp
from tapioca.inputs.svg import SVGInputError, SVGScenario, from_state, reader

from conftest import COARSE, build_scenario, quiet, svg_path


@pytest.fixture(scope="module")
def final_state(tmp_path_factory):
    '''
    Final state of a Mandyoc scenario: the rifting example on the coarse grid (interfaces.txt) with an asthenosphere
    lens in the lithospheric mantle (100-120 km deep, x = 600-1000 km) in lithology.nc and a temperature.nc with
    two times (0 and 5 Myr) whose field changes along x. Return (folder, scenario of the example).
    '''
    scenario = build_scenario("rifting")
    folder = tmp_path_factory.mktemp("preliminary")
    with quiet():
        scenario.write_all(str(folder), plots=False, scripts=False, copy_svg=False)
    x, z = scenario.x, -scenario.z[::-1]  #z from -Lz to 0
    depth = -z
    codes = np.zeros((len(z), len(x)), dtype=np.int8)
    for i, layer in enumerate(scenario.layers):  #bottom to top: each layer covers the ones below it
        codes[depth[:, None] < layer.interface[None, :]] = i
    lens = (depth[:, None] >= 100e3) & (depth[:, None] < 120e3) & (x[None, :] >= 600e3) & (x[None, :] <= 1000e3)
    codes[lens] = 0
    times = [0.0, 5.0]
    xr.Dataset({"lithology": (("time", "z", "x"), np.stack([codes] * 2))},
               coords={"time": times, "z": z, "x": x}).to_netcdf(folder / "lithology.nc")
    T = np.clip(1350.0 * (depth[:, None] - 40e3) / 155e3, 0, None) + 1.0e-4 * x[None, :]
    xr.Dataset({"temperature": (("time", "z", "x"), np.stack([T, T + 10.0]))},
               coords={"time": times, "z": z, "x": x}).to_netcdf(folder / "temperature.nc")
    return folder, scenario


def _svg_from_state(final_state, tmp_path, **kwargs):
    folder, _ = final_state
    svg = tmp_path / "new" / "final_state.svg"
    with quiet():
        agreement = from_state.state_to_svg(str(folder), str(svg), cell=5.0, dx=10.0, dz=10.0, **kwargs)
    return svg, agreement


# Units ##############################################################################

def test_unit_labels_with_weakening_seed(final_state):
    folder, scenario = final_state
    params, _ = tp.read_interfaces(str(folder / "interfaces.txt"))
    assert from_state.unit_labels(params) == [layer.layer_label for layer in scenario.layers]


def test_unit_labels_with_the_seed_found_by_its_C():
    '''Old interfaces.txt (7 parameters, as in the continental rift of Mandyoc): the seed is the mantle with C < 1.'''
    mantle, crust = [3354.0, 9e-12, 2.4168e-15, 3.5, 540000.0, 2.5e-05], [2800.0, 2.86e-10, 8.574e-28, 4.0, 222000.0, 0.0]
    rows = [[1.0, 3378.0, 0.0, 1.393e-14, 3.0, 429000.0, 1.5e-05], [1.0] + mantle, [0.1] + mantle, [1.0] + mantle,
            [1.0] + crust, [1.0, 2700.0, 9.26e-10, 8.574e-28, 4.0, 222000.0, 0.0], [1.0, 1.0, 0.0, 1e-18, 1.0, 0.0, 0.0]]
    params = {key: [row[k] for row in rows] for k, key in enumerate(["C", "rho", "H", "A", "n", "Q", "V"])}
    assert from_state.unit_labels(params) == ["asthenosphere", "lithospheric_mantle", "weak_seed", "lithospheric_mantle",
                                              "lower_crust", "upper_crust", "air"]


# SVG ################################################################################

def test_svg_from_state(final_state, tmp_path):
    svg, agreement = _svg_from_state(final_state, tmp_path, time=5.0)
    assert agreement > 0.99
    setup = reader.read_setup(str(svg))
    assert (setup["Lx"], setup["Lz"]) == (1600e3, 300e3)
    assert setup["experiment"]["preset_temperature_mode"] == "field"
    folder, _ = final_state
    assert (svg.parent / setup["experiment"]["preset_initial_temperature"]).resolve() == folder.resolve()
    assert re.fullmatch(r"5 Myr", setup["experiment"]["preset_time"])


def test_scenario_from_state(final_state, tmp_path):
    '''The SVG of the final state gives its geometry, its LAB (column by column) and its 2D temperature.'''
    folder, original = final_state
    svg, _ = _svg_from_state(final_state, tmp_path, time=5.0)
    with quiet():
        scenario = SVGScenario(str(svg), verbose=False)
        T = scenario.build_temperature()
    x = scenario.x / 1e3
    assert len(reader.layer_indices(scenario.layers, "asthenosphere")) == 2
    #top of the lens: its nodes (100 and 110 km deep) are the cells from 95 to 115 km of the SVG
    np.testing.assert_allclose(scenario.thickness_litho[(x > 620) & (x < 980)], 95e3 - 40e3)
    np.testing.assert_allclose(scenario.thickness_litho[(x < 550) | (x > 1050)], 155e3)  #base of the mantle

    #the temperature of the nodes of temperature.nc at 5 Myr (same grid), with 0 oC in the sticky air
    tb = scenario.temperature_builder
    with xr.open_dataset(folder / "temperature.nc") as dataset:
        expected = dataset.temperature.sel(time=5.0).interp(z=tb.temperature.z, x=tb.temperature.x).values
    air = -tb.temperature.z.values[:, None] < scenario.thickness_sa[None, :]
    np.testing.assert_allclose(T[~air], expected[~air])
    assert np.all(T[air] == 0.0)
    assert scenario.builder.DTree.fields.temperature.attrs["t_potential"] == scenario.TP


def test_field_mode_needs_the_same_domain(final_state, tmp_path):
    folder, _ = final_state
    svg = tmp_path / "other_domain.svg"
    svg.write_text(svg_path("rifting").read_text().replace("Lx = 1600 km", "Lx = 1500 km"))
    with quiet():
        scenario = SVGScenario(str(svg), verbose=False, preset_initial_temperature=str(folder),
                               preset_temperature_mode="field", **COARSE)
    with pytest.raises(SVGInputError, match="differs from the one of the preliminary scenario"):
        scenario.build_temperature()


@pytest.mark.parametrize("settings, message", [({"preset_temperature_mode": "field"}, "needs preset_initial_temperature"),
                                               ({"preset_time": 5}, "only with preset_temperature_mode = field")])
def test_preset_settings_errors(final_state, settings, message):
    folder, _ = final_state
    if("preset_time" in settings):
        settings = {**settings, "preset_initial_temperature": str(folder)}
    with pytest.raises(SVGInputError, match=message):
        build_scenario("rifting", **settings)


def test_cli_tapioca_svg_from_state(final_state, tmp_path, monkeypatch):
    folder, _ = final_state
    svg = tmp_path / "state.svg"
    with quiet():
        from_state.main([str(folder), "-o", str(svg), "--cell", "5", "--dx", "10", "--dz", "10",
                         "--setup", "time_max = 20 Myr"])
    assert reader.read_setup(str(svg))["experiment"]["time_max"] == "20 Myr"
