'''
SVG of a scenario that starts from the final state of a Mandyoc scenario (tapioca-svg-from-state).

The lithology of the final state (lithology.nc of the scenario, last time by default) is rasterized on a regular
grid and vectorized: each unit becomes one SVG path (with its holes), whose description gives its material
parameters, taken from the interfaces.txt of the scenario. The setup shape gives the domain, the resolution and
the experiment, with preset_initial_temperature = <folder of the scenario> and preset_temperature_mode = field,
so tapioca-svg builds the new scenario with the 2D temperature field of the same state:

    tapioca-svg-from-state path/to/scenario -o final_state.svg
    tapioca-svg final_state.svg -o inputs

Units: one per lithology code of interfaces.txt (0 = the bottom layer), labelled from their parameters (see
unit_labels): the bottom one is the asthenosphere and the top one is the sticky air; the weak seed has
weakening_seed = 1 (or the parameters of another code with a smaller C); the other olivine layers are lithospheric
mantle and the other layers are crust (lower_crust and upper_crust when there are two). --labels gives them
instead. Codes with the same label and the same parameters are one unit.

The geometry is that of the final state, so the stack of layers of tapioca may have more slots than the
original scenario (e.g. the asthenosphere that rose in the rift axis), and the lithosphere of each column goes
down to the first asthenosphere below the surface (SVGScenario._find_lab). The accumulated strain is not kept.
'''

import argparse
import os
import sys

import numpy as np
import xarray as xr
from scipy import ndimage
import shapely
from shapely.geometry import box
from shapely.ops import unary_union

from ..._aux_functions import read_interfaces
from ..._variables import DEFAULT_MATERIALS
from .reader import SVGInputError

__all__ = ["unit_labels", "state_to_svg"]

#colors of the units (as in the example SVGs and in the videos of the post-processing); others are grey
COLORS = {'asthenosphere': '#cfe2cd', 'lithospheric_mantle': '#9bc29b', 'weak_seed': '#000000',
          'lower_crust': '#f0d1bc', 'middle_crust': '#eab79f', 'upper_crust': '#e49c7c', 'crust': '#e49c7c',
          'air': '#ffffff'}
GREY = '#b0b0b0'

###############################################################################
# Units
###############################################################################

def rheology_name(params, code):
    '''Material of DEFAULT_MATERIALS with the flow law (A, n, Q) of a code (e.g. 'wet olivine'), or None.'''
    A, n, Q = (float(params[key][code]) for key in ('A', 'n', 'Q'))
    for name, material in DEFAULT_MATERIALS.items():
        if(material['A'] and np.isclose(material['A'], A, rtol=0.01) and material['n'] == n and material['Q'] == Q):
            return name.lower().replace('_', ' ')
    return None

def is_olivine(params, code):
    name = rheology_name(params, code)
    return ('olivine' in name) if(name is not None) else (float(params['rho'][code]) > 3000.0)

def unit_labels(params):
    '''
    Labels of the lithology codes of interfaces.txt (params: {parameter: values of each code, from the bottom}):
    asthenosphere (code 0), air (the top code), weak_seed (weakening_seed = 1, or the parameters of another code
    except a smaller C), lithospheric_mantle (olivine) and crust (lower_crust, middle_crust and upper_crust, from
    the bottom, when there are two or three crustal codes). A label of codes with different parameters gets a
    suffix (_2, _3, ...) from the second one on.
    '''
    ncodes = len(next(iter(params.values())))
    values = np.array([[float(v[code]) for key, v in params.items() if key != 'weakening_seed'] for code in range(ncodes)])
    keys = [key for key in params if key != 'weakening_seed']
    iC = keys.index('C') if('C' in keys) else None

    def is_seed(code):
        if('weakening_seed' in params):
            return float(params['weakening_seed'][code]) == 1
        if(iC is None):
            return False
        others = np.delete(np.arange(len(keys)), iC)
        return any(np.allclose(values[code, others], values[k, others]) and values[code, iC] < values[k, iC]
                   for k in range(1, ncodes - 1) if k != code)

    labels = [None] * ncodes
    labels[0], labels[-1] = 'asthenosphere', 'air'
    crust = []
    for code in range(1, ncodes - 1):
        if(is_seed(code)):
            labels[code] = 'weak_seed'
        elif(is_olivine(params, code)):
            labels[code] = 'lithospheric_mantle'
        else:
            crust.append(code)
    names = {1: ['crust'], 2: ['lower_crust', 'upper_crust'], 3: ['lower_crust', 'middle_crust', 'upper_crust']}
    for k, code in enumerate(crust):
        labels[code] = names[len(crust)][k] if(len(crust) in names) else f'crust_{k + 1}'

    #a label of codes with different parameters: suffix from the second set of parameters on
    for label in set(labels):
        codes = [code for code in range(ncodes) if labels[code] == label]
        distinct = []
        for code in codes:
            match = next((k for k, ref in enumerate(distinct) if np.allclose(values[code], values[ref])), None)
            if(match is None):
                distinct.append(code)
                match = len(distinct) - 1
            if(match > 0):
                labels[code] = f'{label}_{match + 1}'
    return labels

def unit_description(params, code, air=False):
    '''Description of a unit (first line: rheology, then key = value) from the parameters of interfaces.txt.'''
    name = 'air' if(air) else rheology_name(params, code)
    if(name is None):
        name = 'wet olivine' if(is_olivine(params, code)) else 'wet quartz'
        print(f"Warning: flow law of code {code} not in DEFAULT_MATERIALS; '{name}' with the values of interfaces.txt.")
    return '\n'.join([name] + [f'{key} = {float(values[code]):g}' for key, values in params.items()])

###############################################################################
# Raster and SVG
###############################################################################

def remove_small_bodies(grid, min_cells, top_unit=None):
    '''
    Connected regions of a unit with less than min_cells cells, and regions of top_unit (the sticky air) that do
    not reach the top of the domain, take the unit of the nearest cell outside them. Without this, each small
    region would be a layer of the stack of tapioca. Return the grid and the number of cells changed.
    '''
    small = np.zeros(grid.shape, dtype=bool)
    for k in np.unique(grid):
        labels, n = ndimage.label(grid == k)
        remove = np.bincount(labels.ravel()) < min_cells
        if(k == top_unit):
            remove[np.setdiff1d(np.arange(n + 1), labels[0])] = True
        remove[0] = False
        small |= remove[labels]
    if(small.any()):
        idx = ndimage.distance_transform_edt(small, return_distances=False, return_indices=True)
        grid = grid[tuple(idx)]
    return grid, int(small.sum())

def grid_to_polygons(mask, cell_x, cell_z):
    '''Union of the cells of mask (top row = top of the domain) as shapely geometry in km.'''
    boxes = []
    for i, row in enumerate(mask):
        edges = np.flatnonzero(np.diff(np.concatenate([[False], row, [False]]).astype(int)))
        for j0, j1 in zip(edges[::2], edges[1::2]):
            boxes.append(box(j0 * cell_x, i * cell_z, j1 * cell_x, (i + 1) * cell_z))
    return unary_union(boxes)

def svg_path(geometry, ndigits=3):
    '''SVG path data of a (multi)polygon, with the holes in the opposite direction of the exteriors.'''
    polygons = list(geometry.geoms) if hasattr(geometry, 'geoms') else [geometry]
    parts = []
    for polygon in polygons:
        polygon = shapely.geometry.polygon.orient(polygon.simplify(0), sign=1.0)
        for ring in [polygon.exterior] + list(polygon.interiors):
            xy = np.round(np.asarray(ring.coords)[:-1], ndigits)
            parts.append('M ' + ' L '.join(f'{x:g},{y:g}' for x, y in xy) + ' Z')
    return ' '.join(parts)

def write_svg(fname, grid, units, Lx, Lz, setup_lines, comment=''):
    '''
    SVG of tapioca: one path per unit (painted from the bottom unit to the top one) in the layer "lithologies"
    and the "setup" rectangle in the layer "domain". 1 user unit = 1 km.
    units: list of dict(label, title, color, description), one per value of grid.
    '''
    Lx_km, Lz_km = Lx / 1.0e3, Lz / 1.0e3
    cell_x, cell_z = Lx_km / grid.shape[1], Lz_km / grid.shape[0]
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
           f'<!-- {comment} 1 user unit = 1 km; page = model domain ({Lx_km:g} x {Lz_km:g} km). -->',
           f'<svg width="{Lx_km:g}mm" height="{Lz_km:g}mm" viewBox="0 0 {Lx_km:g} {Lz_km:g}" version="1.1" id="svg1"',
           '   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"',
           '   xmlns="http://www.w3.org/2000/svg">',
           '  <g inkscape:groupmode="layer" id="layer1" inkscape:label="lithologies">']
    for k, unit in enumerate(units):
        mask = grid == k
        if(not mask.any()):
            continue
        out += [f'    <path id="{unit["label"]}" inkscape:label="{unit["label"]}" '
                f'd="{svg_path(grid_to_polygons(mask, cell_x, cell_z))}" '
                f'style="fill:{unit["color"]};fill-rule:nonzero;stroke:none">',
                f'      <title id="title_{unit["label"]}">{unit["title"]}</title>',
                f'      <desc id="desc_{unit["label"]}">{unit["description"]}</desc>',
                '    </path>']
    out += ['  </g>',
            '  <g inkscape:groupmode="layer" id="layer_setup" inkscape:label="domain">',
            f'    <rect id="setup" inkscape:label="setup" x="0" y="0" width="{Lx_km:g}" height="{Lz_km:g}" '
            'style="fill:none;stroke:#ff0000;stroke-width:1">',
            '      <title id="title_setup">Model setup</title>',
            '      <desc id="desc_setup">' + '\n'.join(setup_lines) + '</desc>',
            '    </rect>',
            '  </g>',
            '</svg>']
    with open(fname, 'w') as f:
        f.write('\n'.join(out) + '\n')

def state_to_svg(scenario, output='final_state.svg', time=None, labels=None, cell=1.0, min_area=20.0, dx=1.0, dz=1.0,
                 regime='rifting', velocity='1 cm/yr', setup=(), check=None, verbose=True):
    '''
    SVG of a scenario that starts from the final state of a Mandyoc scenario (see the module docstring).

    Parameters
    ----------
    scenario: str
        Folder of the scenario, with lithology.nc, interfaces.txt and temperature.nc.
    output: str
        SVG file.
    time: float
        Time of the state [Myr], the nearest one of lithology.nc (default: the last one). It is also the
        preset_time of the temperature.
    labels: list of str
        Labels of the lithology codes, from the bottom (default: unit_labels).
    cell: float
        Size of the cells of the raster of the SVG [km].
    min_area: float
        Bodies of a unit smaller than min_area [km2] take the unit around them.
    dx, dz: float
        Resolution of the new scenario [km] (setup shape).
    regime, velocity: str
        Experiment of the new scenario (setup shape).
    setup: list of str
        Other lines of the setup shape (e.g. "time_max = 20 Myr").
    check: str
        Figure comparing the raster of the SVG with lithology.nc.

    Returns
    -------
    Agreement of the raster of the SVG with lithology.nc (fraction of the cells).
    '''
    say = print if(verbose) else (lambda *args, **kwargs: None)
    for fname in ['lithology.nc', 'interfaces.txt']:
        if(not os.path.isfile(os.path.join(scenario, fname))):
            raise SVGInputError(f"{fname} not found in {scenario}.")
    params, _ = read_interfaces(os.path.join(scenario, 'interfaces.txt'))
    ncodes = len(next(iter(params.values())))
    labels = list(labels) if(labels is not None) else unit_labels(params)
    if(len(labels) != ncodes):
        raise SVGInputError(f"{len(labels)} labels for the {ncodes} lithology codes of interfaces.txt.")

    with xr.open_dataset(os.path.join(scenario, 'lithology.nc')) as dataset:
        litho = dataset.lithology
        itime = -1 if(time is None) else int(np.abs(litho.time.values - time).argmin())
        litho = litho.isel(time=itime).transpose('z', 'x').sortby('x').load()
    Lx = float(litho.x[-1] - litho.x[0])
    Lz = float(litho.z.max() - litho.z.min())
    t = float(litho.time)
    ref = litho.sortby('z', ascending=False).values.astype(int) #top row = top of the domain
    #a boundary column only of air (e.g. the last column of the lithology.nc of meshes2netcdf.jl) takes the
    #column next to it, otherwise the SVG would have a column of air down to the bottom
    for edge, inner in ((0, 1), (-1, -2)):
        if(np.all(ref[:, edge] == ncodes - 1)):
            ref[:, edge] = ref[:, inner]
            say(f"Column x = {float(litho.x[edge])/1.0e3:g} km of lithology.nc is only air: the one next to it is used.")

    #raster of the SVG: code at the centre of each cell, then one unit per label
    Nx, Nz = int(round(Lx / 1.0e3 / cell)), int(round(Lz / 1.0e3 / cell))
    rows = np.clip(((np.arange(Nz) + 0.5) / Nz * ref.shape[0]).astype(int), 0, ref.shape[0] - 1)
    cols = np.clip(((np.arange(Nx) + 0.5) / Nx * ref.shape[1]).astype(int), 0, ref.shape[1] - 1)
    codes = ref[rows][:, cols]
    unit_names = list(dict.fromkeys(labels))
    code_unit = np.array([unit_names.index(label) for label in labels])
    reference = code_unit[codes]
    units = []
    for label in unit_names:
        code = labels.index(label)
        units.append(dict(label=label, title=label.replace('_', ' ').capitalize(), color=COLORS.get(label, GREY),
                          description=unit_description(params, code, air=(code == ncodes - 1))))

    grid = reference.copy()
    air = code_unit[ncodes - 1]
    for _ in range(3):
        grid, nsmall = remove_small_bodies(grid, max(int(round(min_area / cell**2)), 1), top_unit=air)
        if(nsmall == 0):
            break
        say(f"{nsmall} cells of small bodies replaced by the units around them")

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    folder = os.path.relpath(os.path.abspath(scenario), os.path.dirname(os.path.abspath(output)))
    setup_lines = [f'Lx = {Lx/1.0e3:g} km', f'Lz = {Lz/1.0e3:g} km', f'dx = {dx:g} km', f'dz = {dz:g} km',
                   f'regime = {regime}', f'velocity = {velocity}',
                   f'preset_initial_temperature = {folder}', 'preset_temperature_mode = field']
    if(time is not None):
        setup_lines.append(f'preset_time = {t:g} Myr')
    setup_lines += list(setup)
    write_svg(output, grid, units, Lx, Lz, setup_lines,
              comment=f'Final state of {os.path.abspath(scenario)} (lithology.nc at {t:g} Myr), by tapioca-svg-from-state.')
    agreement = float(np.mean(grid == reference))
    say(f"{output}: {len(units)} units ({', '.join(unit_names)}), raster {Nx} x {Nz}")
    say(f"Agreement with lithology.nc at {t:g} Myr: {100*agreement:.2f}% of the cells")

    if(check):
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap([u['color'] if(u['color'] != '#ffffff') else '#f4f4f4' for u in units])
        extent = (0, Lx/1.0e3, -Lz/1.0e3, 0)
        fig, axs = plt.subplots(3, 1, figsize=(12, 9), constrained_layout=True, sharex=True)
        for ax, data, title in [(axs[0], reference, f'lithology.nc, t = {t:g} Myr'), (axs[1], grid, 'SVG raster'),
                                (axs[2], (grid != reference).astype(float), 'differences')]:
            diff = (title == 'differences')
            ax.imshow(data, extent=extent, aspect='auto', interpolation='nearest', cmap='Greys' if(diff) else cmap,
                      vmin=0 if(diff) else -0.5, vmax=1 if(diff) else len(units) - 0.5)
            ax.set_title(title)
            ax.set_ylabel('z [km]')
        axs[2].set_xlabel('x [km]')
        fig.savefig(check, dpi=150)
        plt.close(fig)
    return agreement

###############################################################################
# Main
###############################################################################

def main(argv=None):
    '''Command line interface: tapioca-svg-from-state path/to/scenario -o final_state.svg'''
    parser = argparse.ArgumentParser(prog='tapioca-svg-from-state', description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('scenario', help='folder of the scenario (lithology.nc, interfaces.txt, temperature.nc)')
    parser.add_argument('-o', '--output', default='final_state.svg', help='SVG file (default: final_state.svg)')
    parser.add_argument('--time', type=float, help='time of the state [Myr] (default: the last one)')
    parser.add_argument('--labels', help='labels of the lithology codes, from the bottom, separated by commas '
                                         '(default: from their parameters)')
    parser.add_argument('--cell', type=float, default=1.0, help='size of the cells of the raster of the SVG [km] (default: 1)')
    parser.add_argument('--min-area', type=float, default=20.0, help='smaller bodies are removed [km2] (default: 20)')
    parser.add_argument('--dx', type=float, default=1.0, help='dx of the new scenario [km] (default: 1)')
    parser.add_argument('--dz', type=float, default=1.0, help='dz of the new scenario [km] (default: 1)')
    parser.add_argument('--regime', default='rifting', help='regime of the new scenario (default: rifting)')
    parser.add_argument('--velocity', default='1 cm/yr', help='velocity of the new scenario (default: 1 cm/yr)')
    parser.add_argument('--setup', action='append', default=[], help='other line of the setup, e.g. "time_max = 20 Myr"')
    parser.add_argument('--check', help='figure comparing the raster of the SVG with lithology.nc')
    args = parser.parse_args(argv)
    try:
        state_to_svg(args.scenario, args.output, time=args.time, labels=args.labels.split(',') if(args.labels) else None,
                     cell=args.cell, min_area=args.min_area, dx=args.dx, dz=args.dz, regime=args.regime,
                     velocity=args.velocity, setup=args.setup, check=args.check)
    except SVGInputError as err:
        sys.exit(f"Error: {err}")


if __name__ == '__main__':
    main()
