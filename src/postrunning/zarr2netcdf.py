#!/usr/bin/env python
'''
Convert the Zarr outputs of Mandyoc (``zarr_output = True``, one ``mandyoc_<step>.zarr`` per printed step, written by
src/zarr_output.cpp of the branch ``zarr_outputs`` of jopmas/mandyoc_jpm) into the same NetCDF4 files that
``meshes2netcdf.jl`` and ``particles2netcdf.jl`` create from the text outputs:

* ``<variable>.nc`` for each grid variable, ``velocity.nc`` with ``vx`` and ``vz``, ``surface.nc``
  (``sp_surface_tracking`` or ``sp_surface_processes``), ``lithology.nc`` (``export_lithology``) and ``times.csv``;
* ``particles_trajectories.nc`` with the particles of the reference step (``print_step_files``).

Same names, dimensions, data types, attributes and compression as the Julia scripts, which keep converting the text
outputs: this script only reads the ``.zarr`` stores and changes nothing for the text outputs. The values come with
full precision from the ``.zarr`` stores instead of the 6 significant digits of the text files.

Usage: python zarr2netcdf.py SCENARIO_DIR [REFERENCE_TIME_STEP]
       (REFERENCE_TIME_STEP of the particle trajectories, default 0, as in particles2netcdf.jl)
'''
import argparse
import csv
import glob
import os
import re
from concurrent.futures import ThreadPoolExecutor

import h5netcdf
import numpy as np
import zarr

try:
    from numba import njit
except ImportError:  # pragma: no cover - numba is in the tapioca environment
    njit = lambda f: f

# ======= the constants of meshes2netcdf.jl and particles2netcdf.jl =======

VARIABLES = ['density', 'viscosity', 'pressure', 'strain', 'strain_rate', 'temperature', 'velocity', 'heat']
AIR_DENSITY_THRESHOLD = -999  # negative: raw data
LITHOLOGY_DATATYPE = np.int8
DEFLATE_LEVEL = 7
PARTICLES_DEFLATE_LEVEL = 5
SELECT_LITHOLOGY = range(0, 128)
THREADS = int(os.environ.get('ZARR2NC_THREADS', os.cpu_count() or 1))

UNITS = {
    'x': 'm', 'y': 'm', 'z': 'm', 'time': 'Myr',
    'density': 'kg/m3', 'viscosity': 'Pa.s', 'pressure': 'Pa', 'strain': 'dimensionless', 'strain_rate': '1/s',
    'temperature': 'deg C', 'velocity': 'm/s', 'surface': 'm', 'heat': 'W/kg', 'thermal_diffusivity': 'm2/s',
    'Phi': 'dimensionless', 'dPhi': '1/s', 'X_depletion': 'dimensionless', 'lithology': 'ID',
}

DTYPES = {
    'x': np.float32, 'y': np.float32, 'z': np.float32, 'time': np.float32,
    'density': np.float32, 'viscosity': np.float64, 'pressure': np.float32, 'strain': np.float64,
    'strain_rate': np.float64, 'temperature': np.float32, 'velocity': np.float64, 'surface': np.float32,
    'heat': np.float64, 'thermal_diffusivity': np.float64, 'Phi': np.float64, 'dPhi': np.float64,
    'X_depletion': np.float64, 'lithology': LITHOLOGY_DATATYPE,
}


def read_param(fpath='param.txt'):
    '''param.txt as a dict (lower-case keys), as read_param of the Julia scripts.'''
    params = {}
    for line in open(fpath):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key_value = line.split('#')[0].replace(' ', '').split('=')
        if len(key_value) == 2:
            params[key_value[0].lower()] = key_value[1]
    return params


def get_all_steps():
    '''Printed steps (mandyoc_<step>.zarr, step >= 0, like the time_<step>.txt files).'''
    steps = []
    for path in glob.glob('mandyoc_*.zarr'):
        match = re.fullmatch(r'mandyoc_(-?\d+)\.zarr', os.path.basename(path))
        if match and int(match.group(1)) >= 0:
            steps.append(int(match.group(1)))
    return np.array(sorted(steps), dtype=np.int32)


def store(step):
    return zarr.open_group(f'mandyoc_{step}.zarr', mode='r')


def read_time(step):
    '''Time [Myr] of the step.'''
    return float(store(step).attrs['time']) / 1e6


def compression(level, shuffle=True):
    return {'compression': 'gzip', 'compression_opts': level, 'shuffle': shuffle}


# ======= grid variables (meshes2netcdf.jl) =======

def read_field(step, variable, dtype):
    '''Field of the step as (z, x) or, for the velocity, (vx, vz); values below 1e-50 set to zero.'''
    data = store(step)[variable][...]
    data[np.abs(data) < 1e-50] = 0
    if variable == 'velocity':
        return data[..., 0].astype(dtype), data[..., 1].astype(dtype)
    return data.astype(dtype)


def read_surface(step):
    return store(step)['surface']['coordinates'][...].astype(DTYPES['surface'])


def create_nc(variable, steps, times, Nx, Nz, Lx, Lz, Nxs=None):
    '''Empty NetCDF file of the variable, with the dimensions, coordinates and attributes of create_nc (Julia).'''
    f = h5netcdf.File(f'{variable}.nc', 'w')
    if variable == 'lithology':
        Nx, Nz = (Nx - 1) * 5 + 1, (Nz - 1) * 5 + 1
    f.dimensions['time'] = len(steps)
    v = f.create_variable('time', ('time',), np.float64, data=np.asarray(times, dtype=np.float64),
                          **compression(DEFLATE_LEVEL))
    v.attrs.update({'units': UNITS['time'], 'long_name': 'time', 'axis': 'T'})
    v = f.create_variable('step', ('time',), np.int32, data=steps, **compression(DEFLATE_LEVEL))
    v.attrs.update({'units': '', 'long_name': 'step'})
    if variable == 'surface':
        f.dimensions['x'] = Nxs
        v = f.create_variable('x', ('x',), DTYPES['x'], data=np.linspace(0, Lx, Nxs).astype(DTYPES['x']),
                              **compression(DEFLATE_LEVEL))
        v.attrs.update({'units': UNITS['x'], 'long_name': 'x', 'axis': 'X'})
        v = f.create_variable(variable, ('time', 'x'), DTYPES[variable], chunks=(1, Nxs), **compression(DEFLATE_LEVEL))
        v.attrs.update({'long_name': variable, 'units': UNITS.get(variable, '-')})
        return f
    f.dimensions['x'] = Nx
    f.dimensions['z'] = Nz
    v = f.create_variable('x', ('x',), DTYPES['x'], data=np.linspace(0, Lx, Nx).astype(DTYPES['x']),
                          **compression(DEFLATE_LEVEL))
    v.attrs.update({'units': UNITS['x'], 'long_name': 'x', 'axis': 'X'})
    v = f.create_variable('z', ('z',), DTYPES['z'], data=np.linspace(-Lz, 0, Nz).astype(DTYPES['z']),
                          **compression(DEFLATE_LEVEL))
    v.attrs.update({'units': UNITS['z'], 'long_name': 'z', 'axis': 'Z'})
    names = ['vx', 'vz'] if variable == 'velocity' else [variable]
    for name in names:
        v = f.create_variable(name, ('time', 'z', 'x'), DTYPES[variable], chunks=(1, Nz, Nx),
                              **compression(DEFLATE_LEVEL))
        if variable == 'velocity':
            v.attrs.update({'units': UNITS[variable], 'long_name': name})
        else:
            v.attrs.update({'long_name': variable, 'units': UNITS.get(variable, '-')})
    return f


def convert_variable(variable, steps, times, Nx, Nz, Lx, Lz):
    '''<variable>.nc from the .zarr stores, one step at a time (the steps are read in parallel threads).'''
    dtype = DTYPES[variable]
    if variable == 'surface':
        f = create_nc(variable, steps, times, Nx, Nz, Lx, Lz, Nxs=read_surface(steps[0]).shape[0])
        reader = lambda step: read_surface(step)[:, 1]
    else:
        f = create_nc(variable, steps, times, Nx, Nz, Lx, Lz)

        def reader(step):
            data = read_field(step, variable, dtype)
            if AIR_DENSITY_THRESHOLD > 0:
                air = read_field(step, 'density', DTYPES['density']) < AIR_DENSITY_THRESHOLD
                for d in (data if variable == 'velocity' else (data,)):
                    d[air] = 0
            return data
    with f, ThreadPoolExecutor(THREADS) as pool:
        for i, data in enumerate(pool.map(reader, steps)):
            if variable == 'velocity':
                f['vx'][i] = data[0]
                f['vz'][i] = data[1]
            else:
                f[variable][i] = data
    print(f'Saved to {variable}.nc\n---------------')


# ======= lithology (convert_litho_to_nc of meshes2netcdf.jl) =======

@njit
def replace_negatives_with_neighbors(mat):
    '''In place, in the order of the Julia function: negative cells (column-major order) take the first
    non-negative of the neighbours (i-1, j), (i+1, j), (i, j-1), (i, j+1), seen as they are updated.'''
    rows, cols = mat.shape
    negatives = []
    for j in range(cols):
        for i in range(rows):
            if mat[i, j] < 0:
                negatives.append((i, j))
    for (i, j) in negatives:
        for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            ni, nj = i + di, j + dj
            if 0 <= ni < rows and 0 <= nj < cols and mat[ni, nj] >= 0:
                mat[i, j] = mat[ni, nj]
                break
    return mat


def read_lithology(step, Nxl, Nzl):
    group = store(step)['lithology']
    ix, iz, layer = group['ix'][...], group['iz'][...], group['layer'][...]
    litho = np.full((Nzl, Nxl), -1, dtype=LITHOLOGY_DATATYPE)
    litho[iz, ix] = layer.astype(LITHOLOGY_DATATYPE)  # later particles overwrite, as in the Julia loop
    replace_negatives_with_neighbors(litho)
    return litho[::-1, :]  # z from the base to the top


def convert_lithology(steps, times, Nx, Nz, Lx, Lz):
    Nxl, Nzl = (Nx - 1) * 5 + 1, (Nz - 1) * 5 + 1
    f = create_nc('lithology', steps, times, Nx, Nz, Lx, Lz)
    with f, ThreadPoolExecutor(THREADS) as pool:
        for i, data in enumerate(pool.map(lambda s: read_lithology(s, Nxl, Nzl), steps)):
            f['lithology'][i] = data
    print('Lithology was converted to NetCDF4')


# ======= particles (particles2netcdf.jl) =======

def read_particles(step):
    group = store(step)['particles']
    coor = group['coordinates'][...]
    return (coor[:, 0].astype(np.float32), coor[:, 1].astype(np.float32), group['itag'][...].astype(np.int32),
            group['layer'][...].astype(LITHOLOGY_DATATYPE))


def convert_particles(steps, times, Lx, Lz, ref_step):
    x, z, ids, litho = read_particles(ref_step)
    keep = (x >= 0) & (x <= Lx) & (z >= -Lz) & (z <= 0) & np.isin(litho, list(SELECT_LITHOLOGY))
    target_ids = np.sort(ids[keep])
    print(f'Total of tracked particles: {len(target_ids)}')
    id_map = np.zeros(target_ids.max() + 1 if len(target_ids) else 1, dtype=np.int64)
    id_map[target_ids] = np.arange(1, len(target_ids) + 1)  # 0: not tracked

    fname = 'particles_trajectories.nc'
    opts = compression(PARTICLES_DEFLATE_LEVEL, shuffle=False)
    with h5netcdf.File(fname, 'w') as f, ThreadPoolExecutor(THREADS) as pool:
        f.dimensions['id'] = len(target_ids)
        f.dimensions['time'] = len(steps)
        v = f.create_variable('id', ('id',), np.int32, data=target_ids, **opts)
        v.attrs.update({'long_name': 'particle_id', 'cf_role': 'trajectory_id'})
        v = f.create_variable('time', ('time',), np.float64, data=np.asarray(times, dtype=np.float64), **opts)
        v.attrs.update({'units': 'myr', 'long_name': 'time', 'standard_name': 'time', 'axis': 'T'})
        f.create_variable('step', ('time',), np.int32, data=steps, **opts)
        for name, std in (('x', 'projection_x_coordinate'), ('z', 'depth')):
            v = f.create_variable(name, ('time', 'id'), np.float32, chunks=(1, max(len(target_ids), 1)), **opts)
            v.attrs.update({'units': 'm', 'long_name': name, 'standard_name': std, 'axis': name.upper()})
        v = f.create_variable('lithology', ('time', 'id'), LITHOLOGY_DATATYPE, chunks=(1, max(len(target_ids), 1)),
                              **opts)
        v.attrs.update({'long_name': 'lithology', 'coordinates': 'time z x'})
        f.attrs.update({'Conventions': 'CF-1.8', 'featureType': 'trajectory', 'description': 'particle trajectories',
                        'reference_timestep': f'{ref_step}'})

        for i, (px, pz, pid, pl) in enumerate(pool.map(read_particles, steps)):
            bx = np.full(len(target_ids), -999.0, dtype=np.float32)
            bz = np.full(len(target_ids), -999.0, dtype=np.float32)
            bl = np.full(len(target_ids), -1, dtype=LITHOLOGY_DATATYPE)
            ok = (pid >= 0) & (pid < len(id_map))
            idx = np.zeros(len(pid), dtype=np.int64)
            idx[ok] = id_map[pid[ok]]
            ok &= idx > 0
            bx[idx[ok] - 1], bz[idx[ok] - 1], bl[idx[ok] - 1] = px[ok], pz[ok], pl[ok]
            f['x'][i], f['z'][i], f['lithology'][i] = bx, bz, bl
    print(f'Saved to {fname}')


# ======= main =======

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('scenario_dir')
    parser.add_argument('reference_step', nargs='?', type=int, default=0,
                        help='reference time step of the particle trajectories (default: 0)')
    args = parser.parse_args(argv)
    os.chdir(args.scenario_dir)

    params = read_param('param.txt')
    dims = int(params.get('dimensions', '2'))
    if dims != 2:
        raise SystemExit('zarr2netcdf.py converts 2D scenarios only (as meshes2netcdf.jl).')
    if params.get('iterative', 'direct') == 'iterative' or params.get('nondimensionalization', 'False') == 'True':
        for key in DTYPES:
            if key != 'lithology':
                DTYPES[key] = np.float64
        print('You did run a non dimensional model, all datatypes changed to Float64.')

    Nx, Nz = int(params['nx']), int(params['nz'])
    Lx, Lz = float(params['lx']), float(params['lz'])

    steps = get_all_steps()
    if len(steps) == 0:
        raise SystemExit('No mandyoc_<step>.zarr found: run meshes2netcdf.jl for the text outputs.')
    times = np.array([read_time(s) for s in steps], dtype=DTYPES['time'])
    with open('times.csv', 'w', newline='') as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(['step', 'time_myr'])
        writer.writerows([[int(s), str(t)] for s, t in zip(steps, times)])
    print(f'{len(steps)} time steps were found, from {times[0]} Myr [{steps[0]}] up to {times[-1]} Myr [{steps[-1]}].')

    variables = list(VARIABLES)
    if params.get('sp_surface_tracking', 'False') == 'True' or params.get('sp_surface_processes', 'False') == 'True':
        variables.append('surface')
        print('Surface was tracked.')
    if params.get('magmatism', 'off') == 'on':
        variables += ['Phi', 'dPhi', 'X_depletion']
        print('Magmatism was on.')
    if params.get('export_thermal_diffusivity', 'False') == 'True':
        variables.append('thermal_diffusivity')
        print('Thermal diffusivity (kappa) was exported.')

    for variable in dict.fromkeys(variables):
        print(f'Converting: {variable}')
        convert_variable(variable, steps, times, Nx, Nz, Lx, Lz)

    if params.get('export_lithology', 'False') == 'True':
        print('Lithology grid was exported.')
        convert_lithology(steps, times, Nx, Nz, Lx, Lz)

    if params.get('print_step_files', 'True') == 'True' and 'particles' in store(steps[0]):
        if args.reference_step not in steps:
            raise SystemExit(f'Reference step {args.reference_step} has no mandyoc_{args.reference_step}.zarr.')
        convert_particles(steps, times, Lx, Lz, args.reference_step)

    print('Finished')
    print(f'Variables converted: {"; ".join(dict.fromkeys(variables))}')
    print(f'Compression level: {DEFLATE_LEVEL}')


if __name__ == '__main__':
    main()
