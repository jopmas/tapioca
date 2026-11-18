'''
Velocity field of the punch of Mandyoc (examples/punch/create_input_files_punch.ipynb), for custom_veloc_field of
punch.svg: a rigid punch moves down with vz = -1 (non-dimensional, as in the notebook) at the top of the domain, for
|x - Lx/2| <= HALF_WIDTH; the velocity is zero everywhere else.

write_files writes bcv_0.txt, the boundary conditions of each node read by Mandyoc with bcv_extern = True
(1: free, 0: fixed): the velocity is fixed under the punch, vx is fixed at the sides and at the base and vz at the
base, and the other nodes are free.
'''

import os

import numpy as np

from tapioca.inputs import VelocityFieldBuilder

HALF_WIDTH = 0.08  # half width of the punch [m]
VZ_PUNCH = -1.0    # velocity of the punch [m/s], non-dimensional in the notebook


def _under_punch(scenario, x):
    return np.abs(x - scenario.Lx / 2) <= HALF_WIDTH + 1.0e-9 * scenario.Lx


def build_velocity(scenario):
    '''Velocity of the punch at the top of the domain (VelocityFieldBuilder in m/s).'''
    v = VelocityFieldBuilder(scenario.builder, units='m/s')
    tol = 1.0e-9 * scenario.Lx
    v.set_region('top', [scenario.Lx / 2 - HALF_WIDTH - tol, scenario.Lx / 2 + HALF_WIDTH + tol], vx=0.0, vz=VZ_PUNCH)
    return v


def write_files(scenario, output_dir):
    '''bcv_0.txt of the notebook (bcv_extern = True); return the names of the files written.'''
    X, Z = scenario.X, scenario.Z  # Z: depth, from the bottom (Lz) to the top (0), as in the notebook
    vx, vz = np.ones_like(X), np.ones_like(X)
    top = np.isclose(Z, 0.0) & _under_punch(scenario, X)
    vx[top], vz[top] = 0.0, 0.0
    vx[np.isclose(X, 0.0) | np.isclose(X, scenario.Lx) | np.isclose(Z, scenario.Lz)] = 0.0
    vz[np.isclose(Z, scenario.Lz)] = 0.0
    bcv = np.column_stack([vx.ravel(), vz.ravel()]).ravel()
    np.savetxt(os.path.join(output_dir, 'bcv_0.txt'), bcv, header="v1\nv2\nv3\nv4", fmt="%.4f")
    return ['bcv_0.txt']


def infos(scenario):
    '''Lines of infos_<folder>.txt about the velocity field.'''
    return [f"punch: vz = {VZ_PUNCH:g} at the top for |x - Lx/2| <= {HALF_WIDTH:g} m; bcv_0.txt with bcv_extern = True"]
