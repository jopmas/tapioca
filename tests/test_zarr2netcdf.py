'''zarr2netcdf.py (src/postrunning): .nc files from the Zarr outputs of Mandyoc (zarr_output = True).'''
import importlib.util
import os

import numpy as np
import pytest
import xarray as xr

zarr = pytest.importorskip("zarr")

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "src", "postrunning", "zarr2netcdf.py")
NX, NZ, LX, LZ = 5, 4, 4000.0, 3000.0


def load_script():
    spec = importlib.util.spec_from_file_location("zarr2netcdf", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_step(folder, step, time, shift):
    '''A mandyoc_<step>.zarr store as written by src/zarr_output.cpp of mandyoc_jpm (Zarr v2).'''
    root = zarr.open_group(str(folder / f"mandyoc_{step}.zarr"), mode="w", zarr_format=2)
    root.attrs["time"] = time
    field = np.arange(NZ * NX, dtype=np.float64).reshape(NZ, NX) + shift
    for name in ["density", "viscosity", "pressure", "strain", "strain_rate", "temperature", "heat"]:
        root.create_array(name, data=field, fill_value=None)
    root.create_array("velocity", data=np.stack([field, -field], axis=-1), fill_value=None)
    surface = root.create_group("surface")
    surface.create_array("coordinates", data=np.column_stack([np.linspace(0, LX, 7), np.full(7, shift)]),
                         fill_value=None)
    # one particle per node of the lithology grid (5x finer), but none at (ix, iz) = (0, 0)
    ix, iz = np.meshgrid(np.arange((NX - 1) * 5 + 1), np.arange((NZ - 1) * 5 + 1))
    ix, iz = ix.ravel()[1:], iz.ravel()[1:]
    litho = root.create_group("lithology")
    litho.create_array("ix", data=ix.astype(np.int32), fill_value=None)
    litho.create_array("iz", data=iz.astype(np.int32), fill_value=None)
    litho.create_array("layer", data=(iz > 7).astype(np.int32), fill_value=None)
    particles = root.create_group("particles")
    ids = np.array([10003, 10001, 10002] if step else [10001, 10002, 10003], dtype=np.int32)
    particles.create_array("coordinates", data=np.column_stack([ids - 10000.0 + shift, -ids + 10000.0]),
                           fill_value=None)
    particles.create_array("itag", data=ids, fill_value=None)
    particles.create_array("layer", data=np.ones(3, dtype=np.int32), fill_value=None)


@pytest.fixture
def converted(tmp_path):
    (tmp_path / "param.txt").write_text(
        f"nx = {NX}\nnz = {NZ}\nlx = {LX}\nlz = {LZ}\nsp_surface_tracking = True\nexport_lithology = True\n"
        "print_step_files = True\n")
    write_step(tmp_path, 0, 0.0, 0.0)
    write_step(tmp_path, 50, 2.5e5, 1.0)
    cwd = os.getcwd()
    try:
        load_script().main([str(tmp_path), "50"])
    finally:
        os.chdir(cwd)
    return tmp_path


def test_grid_variables(converted):
    with xr.open_dataset(converted / "temperature.nc") as ds:
        assert ds.temperature.dims == ("time", "z", "x") and ds.temperature.dtype == np.float32
        assert list(ds.step.values) == [0, 50] and np.allclose(ds.time, [0, 0.25])
        assert np.allclose(ds.z, np.linspace(-LZ, 0, NZ)) and np.allclose(ds.x, np.linspace(0, LX, NX))
        assert np.allclose(ds.temperature[1], np.arange(NZ * NX).reshape(NZ, NX) + 1)
    with xr.open_dataset(converted / "velocity.nc") as ds:
        assert np.allclose(ds.vz, -ds.vx)
    assert (converted / "times.csv").read_text().splitlines() == ["step,time_myr", "0,0.0", "50,0.25"]


def test_surface_and_lithology(converted):
    with xr.open_dataset(converted / "surface.nc") as ds:
        assert ds.surface.dims == ("time", "x") and np.allclose(ds.surface[1], 1.0)
    with xr.open_dataset(converted / "lithology.nc") as ds:
        litho = ds.lithology.values
        assert litho.shape == (2, (NZ - 1) * 5 + 1, (NX - 1) * 5 + 1) and litho.dtype == np.int8
        # z from the base: iz counts from the top; the empty node (0, 0) takes its neighbour
        assert (litho[:, -8:, :] == 0).all() and (litho[:, :-8, :] == 1).all()


def test_particle_trajectories(converted):
    with xr.open_dataset(converted / "particles_trajectories.nc") as ds:
        assert list(ds.id.values) == [10001, 10002, 10003] and ds.attrs["reference_timestep"] == "50"
        assert np.allclose(ds.x.values, [[1, 2, 3], [2, 3, 4]])  # sorted by id, also when the order changes
