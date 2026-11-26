#!/usr/bin/env python
'''
Convert the Zarr outputs of Mandyoc into the NetCDF4 files that ``meshes2netcdf.jl`` and ``particles2netcdf.jl``
create from the text outputs.

With ``zarr_output = True`` in ``param.txt`` (src/zarr_output.cpp of the branch ``zarr_outputs`` of
jopmas/mandyoc_jpm), Mandyoc writes one Zarr v2 store per variable, in the layout of mandyocjll:
``density.zarr``, ..., ``velocity.zarr`` (``vx``, ``vz``), ``lithology.zarr``, ``particles_trajectories.zarr`` and
``surface.zarr``, with the names, dimensions, units and types of the ``.nc`` files. Each ``<variable>.zarr`` becomes
``<variable>.nc``, with the compression of the Julia scripts (deflate 7 and shuffle; 5 without shuffle for the
trajectories) and ``time`` in float64 as there, and ``times.csv`` is written. The values come with full precision
instead of the 6 significant digits of the text files.

The Julia scripts are unchanged and keep converting the text outputs.

Usage: python zarr2netcdf.py SCENARIO_DIR [variable ...]   (default: every <variable>.zarr of the folder)
'''
import argparse
import csv
import os
from concurrent.futures import ThreadPoolExecutor

import h5netcdf
import numpy as np
import zarr

DEFLATE_LEVEL = 7
PARTICLES_DEFLATE_LEVEL = 5
THREADS = int(os.environ.get('ZARR2NC_THREADS', os.cpu_count() or 1))
STORES = ['density', 'viscosity', 'pressure', 'strain', 'strain_rate', 'temperature', 'velocity', 'heat',
          'thermal_diffusivity', 'Phi', 'dPhi', 'X_depletion', 'surface', 'lithology', 'particles_trajectories']
"""list: Stores of Mandyoc, in the order of meshes2netcdf.jl and particles2netcdf.jl."""


def attributes(array):
    return {k: v for k, v in array.attrs.items() if k != '_ARRAY_DIMENSIONS'}


def convert_store(name):
    '''<name>.zarr -> <name>.nc, one time at a time (read in parallel threads).'''
    group = zarr.open_group(f'{name}.zarr', mode='r')
    arrays = {key: group[key] for key in sorted(group.array_keys())}
    dims = {key: list(a.attrs['_ARRAY_DIMENSIONS']) for key, a in arrays.items()}
    coords = [key for key in arrays if dims[key] == [key]]
    order = {'time': 0, 'step': 1, 'id': 2, 'x': 3, 'z': 4}
    coords.sort(key=lambda k: order.get(k, 9))
    data = [key for key in arrays if key not in coords]
    level, shuffle = (PARTICLES_DEFLATE_LEVEL, False) if name == 'particles_trajectories' else (DEFLATE_LEVEL, True)
    opts = {'compression': 'gzip', 'compression_opts': level, 'shuffle': shuffle}
    nt = arrays['time'].shape[0]

    with h5netcdf.File(f'{name}.nc', 'w') as f:
        f.attrs.update(dict(group.attrs))
        sizes = {}
        for key, a in arrays.items():
            sizes.update(zip(dims[key], a.shape))
        for d, n in sizes.items():
            f.dimensions[d] = n
        for key in coords:
            values = arrays[key][...]
            values = values.astype(np.float64) if key == 'time' else values
            f.create_variable(key, (key,), values.dtype, data=values, **opts).attrs.update(attributes(arrays[key]))
        for key in data:
            a = arrays[key]
            if dims[key][0] != 'time' or len(dims[key]) == 1:  # e.g. step (time)
                f.create_variable(key, tuple(dims[key]), a.dtype, data=a[...], **opts).attrs.update(attributes(a))
                continue
            v = f.create_variable(key, tuple(dims[key]), a.dtype, chunks=(1,) + tuple(a.shape[1:]), **opts)
            v.attrs.update(attributes(a))
            with ThreadPoolExecutor(THREADS) as pool:
                for i, values in enumerate(pool.map(lambda t: a[t], range(nt))):
                    v[i] = values
    print(f'Saved to {name}.nc ({nt} times)')


def write_times_csv(name):
    '''times.csv (step, time_myr), as meshes2netcdf.jl.'''
    group = zarr.open_group(f'{name}.zarr', mode='r')
    with open('times.csv', 'w', newline='') as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(['step', 'time_myr'])
        writer.writerows([[int(s), str(np.float32(t))] for s, t in zip(group['step'][...], group['time'][...])])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('scenario_dir')
    parser.add_argument('variables', nargs='*', help='stores to convert (default: every <variable>.zarr)')
    args = parser.parse_args(argv)
    os.chdir(args.scenario_dir)

    names = args.variables or [n for n in STORES if os.path.isdir(f'{n}.zarr')]
    names += [] if args.variables else sorted(p[:-5] for p in os.listdir('.')
                                              if p.endswith('.zarr') and p[:-5] not in STORES and os.path.isdir(p))
    if not names:
        raise SystemExit('No <variable>.zarr found: run meshes2netcdf.jl for the text outputs.')
    write_times_csv(names[0])
    for name in names:
        print(f'Converting: {name}')
        convert_store(name)
    print('Finished')
    print(f'Variables converted: {"; ".join(names)}')


if __name__ == '__main__':
    main()
