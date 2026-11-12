"""
Shared fixtures of the tests of tapioca.

Run the tests from the root of the repository with:

    pytest
"""

import contextlib
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # figures without a display

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
SVG_EXAMPLES = EXAMPLES / "svg_scenarios"
CONTINENTAL_RIFT = EXAMPLES / "continental_rift"

# example folder -> (regime, moving boundary) written in its SVG
SVG_SCENARIOS = {
    "rifting": ("rifting", "right"),
    "rifting_symmetric": ("rifting", "symmetric"),
    "rifting_with_basin": ("rifting", "right"),
    "compression": ("compression", "left"),
    "accordion": ("accordion", "right"),
    "quiescence": ("quiescence", "right"),
    "rifting_sedimentation": ("rifting", "symmetric"),
}

# coarse grid and short conductive model, so that the SVG scenarios are built in a fraction of a second
COARSE = dict(dx=10, dz=10)
CONDUCTIVE_TIME = 50.0e6  # yr


def svg_path(name):
    '''Path of the SVG of an example in examples/svg_scenarios.'''
    return SVG_EXAMPLES / name / f"{name}.svg"


@contextlib.contextmanager
def quiet():
    '''Silence the messages printed by the builders.'''
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def build_scenario(name, **kwargs):
    '''SVGScenario of an example on the coarse grid, with a short conductive model.'''
    from tapioca.inputs.svg import SVGScenario
    with quiet():
        scenario = SVGScenario(str(svg_path(name)), verbose=False, **{**COARSE, **kwargs})
    scenario.conductive_time = CONDUCTIVE_TIME
    return scenario


@pytest.fixture(scope="session")
def continental_rift():
    '''MandyocScen of examples/continental_rift with all the outputs loaded.'''
    import tapioca as tp
    with quiet():
        return tp.MandyocScen(str(CONTINENTAL_RIFT),
                              variables=["density", "temperature", "velocity", "viscosity", "strain_rate"],
                              load_lithology=True, load_surface=True, load_particles=True)
