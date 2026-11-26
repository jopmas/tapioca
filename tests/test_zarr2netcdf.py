'''zarr2netcdf.py (src/postrunning): .nc files from the Zarr outputs of Mandyoc (zarr_output = True), one
<variable>.zarr store per variable in the layout of mandyocjll.'''
import importlib.util
import os

import h5py
import numpy as np
import pytest
import xarray as xr

zarr = pytest.importorskip("zarr")

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "src", "postrunning", "zarr2netcdf.py")
NX, NZ, LX, LZ, NT = 5, 4, 4000.0, 3000.0, 3


def load_script():
    spec = importlib.util.spec_from_file_location("zarr2netcdf", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def array(group, name, data, dims, chunks=None, **attrs):
    a = group.create_array(name, data=data, chunks=chunks or data.shape, fill_value=None)
    a.attrs.update({"_ARRAY_DIMENSIONS": dims, **attrs})


def write_store(folder, name, variables, extra_dims=(), group_attrs=None):
    '''A <name>.zarr store as written by src/zarr_output.cpp of mandyoc_jpm (Zarr v2).'''
    group = zarr.open_group(str(folder / f"{name}.zarr"), mode="w", zarr_format=2)
    group.attrs.update(group_attrs or {})
    array(group, "time", np.array([0.0, 0.5, 1.0], dtype=np.float32), ["time"], units="Myr", long_name="time", axis="T")
    array(group, "step", np.array([0, 50, 100], dtype=np.int32), ["time"], units="", long_name="step")
    for coord, n, a, b in extra_dims:
        array(group, coord, np.linspace(a, b, n).astype(np.float32), [coord], units="m", long_name=coord,
              axis=coord.upper())
    for var, (data, dims, attrs) in variables.items():
        array(group, var, data, dims, chunks=(1,) + data.shape[1:], **attrs)


@pytest.fixture
def converted(tmp_path):
    field = np.arange(NT * NZ * NX, dtype=np.float32).reshape(NT, NZ, NX)
    mesh = (("x", NX, 0, LX), ("z", NZ, -LZ, 0))
    write_store(tmp_path, "temperature", {"temperature": (field, ["time", "z", "x"], {"units": "deg C"})}, mesh)
    write_store(tmp_path, "velocity", {"vx": (field.astype(np.float64), ["time", "z", "x"], {"units": "m/s"}),
                                       "vz": (-field.astype(np.float64), ["time", "z", "x"], {"units": "m/s"})}, mesh)
    group_attrs = {"Conventions": "CF-1.8", "featureType": "trajectory", "reference_timestep": "0"}
    write_store(tmp_path, "particles_trajectories",
                {"x": (np.ones((NT, 2), dtype=np.float32), ["time", "id"], {"units": "m"}),
                 "lithology": (np.full((NT, 2), -1, dtype=np.int8), ["time", "id"], {})}, group_attrs=group_attrs)
    group = zarr.open_group(str(tmp_path / "particles_trajectories.zarr"), mode="a", zarr_format=2)
    array(group, "id", np.array([10000, 10001], dtype=np.int32), ["id"], long_name="particle_id",
          cf_role="trajectory_id")
    cwd = os.getcwd()
    try:
        load_script().main([str(tmp_path)])
    finally:
        os.chdir(cwd)
    return tmp_path


def test_mesh_store(converted):
    with xr.open_dataset(converted / "temperature.nc") as ds:
        assert ds.temperature.dims == ("time", "z", "x") and ds.temperature.dtype == np.float32
        assert ds.time.dtype == np.float64 and list(ds.step.values) == [0, 50, 100]
        assert np.allclose(ds.z, np.linspace(-LZ, 0, NZ)) and ds.temperature.attrs["units"] == "deg C"
        assert np.allclose(ds.temperature[2], np.arange(NZ * NX).reshape(NZ, NX) + 2 * NZ * NX)
    with xr.open_dataset(converted / "velocity.nc") as ds:
        assert np.allclose(ds.vz, -ds.vx)
    with h5py.File(converted / "temperature.nc") as f:  # compression of meshes2netcdf.jl
        assert (f["temperature"].compression_opts, f["temperature"].shuffle) == (7, True)
    assert (converted / "times.csv").read_text().splitlines() == ["step,time_myr", "0,0.0", "50,0.5", "100,1.0"]


def test_trajectories(converted):
    with xr.open_dataset(converted / "particles_trajectories.nc") as ds:
        assert ds.x.dims == ("time", "id") and list(ds.id.values) == [10000, 10001]
        assert ds.attrs["featureType"] == "trajectory" and ds.attrs["reference_timestep"] == "0"
        assert ds.lithology.dtype == np.int8
    with h5py.File(converted / "particles_trajectories.nc") as f:  # compression of particles2netcdf.jl
        assert (f["x"].compression_opts, f["x"].shuffle) == (5, False)


def test_only_the_given_variables(converted):
    for f in converted.glob("*.nc"):
        f.unlink()
    cwd = os.getcwd()
    try:
        load_script().main([str(converted), "velocity"])
    finally:
        os.chdir(cwd)
    assert sorted(p.name for p in converted.glob("*.nc")) == ["velocity.nc"]
