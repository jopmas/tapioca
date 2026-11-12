'''
Builder of the input files of a Mandyoc scenario drawn in an SVG file.

The geometry and the materials come from the SVG (see tapioca.inputs.svg.reader). The shape
labelled "setup" gives the domain, the resolution and, optionally, the numerical
experiment:

    Lx = 1600 km
    Lz = 400 km
    dx = 2 km
    dz = 2 km
    regime = accordion                         (rifting, compression, accordion, quiescence)
    velocity = 1 cm/yr
    moving_boundary = right                    (left, right or symmetric)
    quiescence_time = 0 Myr
    time_max = 500 Myr
    magmatism = off
    magmatism_extraction = off
    magmatism_Ts0 = 1080
    heat_capacity = 1250                       (J/kg/K)
    thermal_expansion_coefficient = 3.28e-5    (1/K)
    gravity_acceleration = 10                  (m/s2)
    preset_initial_temperature = ../preliminary
    custom_veloc_field = build_veloc.py        (Python script with build_velocity(scenario))
    temperature_solver = 2D                    (2D or columns)
    min_thickness = 1 km
    water_layer = water                        (drawn below the air, or a template outside the domain)
    sp_surface_processes = True                (surface processes, see below)
    sp_mode = theunissen_sedimentation
    sp_format = ggciag                         (ggciag or SP_vc: files of the version of Mandyoc)
    sediment_layer = sediments                 (unit of the deposited layers, usually a template outside the domain)
    ti_sed = 0 Myr
    tf_sed = 70 Myr                            (default: time_max)
    dt_sed = 1 Myr
    sedimentation_rate = 25                    (or in steps: 0 Myr: 25, 20 Myr: 50, 40 Myr: 10)
    Ld = 35 km
    base_level = -500 m
    strain_sed = 0.25                          (and the other parameters of SP_PARAMS, written to param.txt)

Sedimentation: with sp_surface_processes = True and a sp_mode that deposits sediments (SEDIMENTATION_MODES), one
layer of zero thickness with the material of sediment_layer is added at the top of the solid layers for each
interval of dt_sed between ti_sed and tf_sed. sedimentation_rate [m2/yr] may change in time: "t1: rate1, t2: rate2,
..." (times in Myr, or with units) is constant in steps, each rate from its time on. The unit of sediment_layer does
not need to be drawn in the domain: a shape drawn outside the setup shape gives its material (a template), and the
same holds for the water of water_layer (a layer of zero thickness below the sticky air) and for the magmatic_layer
(a layer of zero thickness at the top of the solid layers; without magmatic_layer, the basalt of tapioca). With
base_level, base_level.txt keeps the base level of Mandyoc at that height [m] (variable_base_level = True).

The files depend on the version of Mandyoc (sp_format):
  - ggciag (default), the Mandyoc of ggciag (https://github.com/ggciag/mandyoc): sediment_layers.txt (the time and the
    ID of each layer; written with every sp_mode, since Mandyoc reads it whenever sp_surface_processes is True),
    sedimentation_rate.txt (the time and the rate of each change) and base_level.txt, each with the number of lines
    first. Ld is fixed at 35 km in its code and it has no water layer, so water_interface is not written;
  - SP_vc, the Mandyoc with the sedimentation of Theunissen (mandyoc_tests, of the group): SP_vc.txt (one line per
    change: time [Myr], ID of the active layer, sedimentation_rate and Ld [m]; a change of the rate between the times
    of two layers adds a line with the layer that is active then), base_level.txt without the number of lines and
    water_interface in param.txt. The Mandyoc of ggciag stops at water_interface, which it does not know.
With both, a tf_sed before time_max adds a last change at tf_sed with rate 0, which stops the sedimentation.

Each setting may also be given to SVGScenario (or on the command line, tapioca-svg),
which overrides the SVG; without it, the default is used.

Example
-------
    from tapioca.inputs.svg import SVGScenario
    scenario = SVGScenario("rifting.svg", regime="accordion", quiescence_time=0)
    scenario.write_all("my_scenario")

The scenario is built on the generic builders of tapioca: scenario.builder is a MandyocBuilder with the
interfaces and the material fields of the SVG, and build_temperature and build_velocity use a
TemperatureFieldBuilder and a VelocityFieldBuilder. The input files are exported by scenario.builder,
so interfaces.txt has the 13 parameters of tapioca.DEFAULT_MATERIAL.
'''

import os, glob, shutil, platform, zipfile, importlib.util

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ..._variables import SEC_PER_YEAR, PARAMS_TYPES
from ..._aux_functions import find_nearest, read_preset_temperature, adapt_preset_temperature
from ..generic_input_classes import TemperatureFieldBuilder, VelocityFieldBuilder, VerboseExtension, GenericTimeConditions
from ..checks import SVGChecks
from ._variables import (REGIMES, MANTLE_RHEOLOGIES, SYSTEMS, DEFAULT_PARAMS, SP_MODES, SEDIMENTATION_MODES, SP_PARAMS, SP_FORMATS,
                         ORGANIZE_SCRIPT, RUN_SCRIPT,
                         POST_PROCESSING_SCRIPT, NO_POST_PROCESSING_SCRIPT,
                         RUN_LAUNCHER, PETSC_ENV, PETSC_ENV_LOOKUP, ZIP_SCRIPT, HYPATIA, HYPATIA_SCRIPT,
                         HYPATIA_CONDA_SCRIPT)
from .materials import MandyocLayer, BASALT
from ... import __version__
from .reader import (SVGInputError, SETUP_LABEL, MIN_THICKNESS, EXPERIMENT_KEYS, SVG_VERSION, INKSCAPE_VERSION,
                     svg_versions, read_setup, layers_from_svg, layers_to_builder, points_from_spacing, parse_quantity,
                     layer_indices, slot_thickness, layer_thickness, unit_layer)

# (setting, minimum, value must be greater than the minimum): used in param.txt and in the initial temperature
THERMAL_SETTINGS = (("heat_capacity", 0.0, True),
                    ("thermal_expansion_coefficient", 0.0, False),
                    ("gravity_acceleration", 0.0, True))

__all__ = ["SVGScenario", "REGIMES"]

###############################################################################
# Scenario builder
###############################################################################

class SVGScenario(SVGChecks, VerboseExtension):
    '''
    Mandyoc scenario drawn in an SVG file.

    The geometry, the materials and the domain come from the SVG (see tapioca.inputs.svg.reader). The settings
    of the numerical experiment are taken from the arguments, then from the description of the setup
    shape of the SVG and finally from the defaults below.

    Parameters
    ----------
    svg_file : str
        SVG file with the scenario.
    dx, dz : float, optional
        Horizontal and vertical resolution [km] (default: from the SVG).

    Settings of the numerical experiment (keyword arguments, see reader.EXPERIMENT_KEYS); each one overrides the
    setting of the setup shape of the SVG. Numbers may be given with units, as in the SVG (e.g. '1 cm/yr'):

    regime : str, optional
        'rifting' (extension), 'compression' (convergence), 'accordion' (extension, convergence,
        quiescence and extension again, as in M300) or 'quiescence' (no velocity at the boundaries).
        Default: rifting.
    velocity : float, optional
        Velocity at the moving boundary [cm/yr] (with symmetric, the sum of the velocities of both boundaries).
        Default: 1 cm/yr.
    moving_boundary : str, optional
        'left', 'right' or 'symmetric' (both boundaries move with velocity/2 and the reference frame is at the
        centre of the domain). Default: right.
    quiescence_time : float, optional
        Duration of the quiescence [Myr] in the accordion and quiescence regimes. Default: 0 Myr.
    time_max : float, optional
        Duration of the simulation [Myr]. Default: 60 Myr in rifting and compression, the end of the
        second extension in accordion (which time_max sets when given) and quiescence_time (or 1000 Myr)
        in quiescence.
    magmatism, magmatism_extraction : str, optional
        'on' or 'off'. Default: off. With both on, a basalt layer of zero thickness receives the extracted
        melt (see magmatic_layer).
    magmatism_Ts0 : float, optional
        magmatism_Ts0 of param.txt [oC]. Default: 1080.
    heat_capacity, thermal_expansion_coefficient, gravity_acceleration : float, optional
        Heat capacity [J/kg/K], thermal expansion coefficient [1/K] and gravity (positive) [m/s2]. Each value is
        used in param.txt, in the initial temperature and in the yield strength envelopes, so they always agree.
        Plain numbers, without units. Default: 1250, 3.28e-5 and 10 (tapioca.inputs.svg._variables.DEFAULT_PARAMS).
    water_layer : str or int, optional
        Label (or ID from the bottom) of the unit drawn as water, right below the sticky air. The surface (top of
        the solid layers) is then its base. Default: none (no water).
    magmatic_layer : str or int, optional
        Label (or ID from the bottom) of the layer that receives the melt extracted by Mandyoc (magmatic_layer of
        param.txt). A unit drawn only outside the domain (a template, e.g. a basalt with its own material) becomes a
        layer of zero thickness at the top of the solid layers, below the sediments and the water. Default: with magmatism and magmatism_extraction on, a basalt layer of zero thickness added at
        the top of the solid layers (below the water if there is water, otherwise below the sticky air); otherwise
        the top solid layer. A layer given here is used as it is, so no basalt layer is added.
    preset_initial_temperature : str, optional
        Folder (with temperature.nc, or param.txt and temperature_*.txt) or NetCDF file of a preliminary
        scenario (see tapioca.read_preset_temperature): the horizontal mean of its final
        temperature, adapted to the lithosphere of the SVG, is the initial temperature. A relative path in the
        SVG is relative to the folder of the SVG. Default: none (conductive model).
    custom_veloc_field : str, optional
        Python script that builds the velocity field instead of the kinematic regime (see build_velocity). A
        relative path in the SVG is relative to the folder of the SVG. Default: none (kinematic regime).
    temperature_solver : str, optional
        Steady-state conduction in the lithosphere of the conductive model: '2D' (the whole lithosphere at once,
        so heat also flows sideways where the layers change along x) or 'columns' (each column alone).
        Default: 2D.

    Other parameters:

    asthenosphere_label, seed_label : str
        Labels of the asthenosphere and of the weak seed in the SVG.
    sampling : float
        Vertical sampling used to locate the interfaces in the SVG [m].
    min_thickness : float, optional
        Pieces of a unit thinner than min_thickness [km] in every column, when the unit is thicker
        elsewhere, are removed (slivers left by overlapping shapes; 0 keeps them). Default: from the setup
        shape of the SVG, or 1 km.
    verbose : bool
        Print the steps.

    Attributes
    ----------
    layers : list of MandyocLayer (bottom to top); layer.interface is the depth of the layer base [m]
    x, z : coordinates of the nodes [m] (z is the depth, from the bottom to the top of the domain)
    T : initial temperature (Nz, Nx), after build_temperature()
    VX, VZ : initial velocity (Nz, Nx), after build_velocity()
    builder : MandyocBuilder (tapioca.inputs) with the interfaces and the fields of the scenario; the
              input files are exported by it, so it can be changed by script before write_all
    temperature_builder, velocity_builder : TemperatureFieldBuilder and VelocityFieldBuilder of the scenario,
              after build_temperature() and build_velocity()
    '''

    #Initial temperature (the thermal diffusivity k of each layer comes from its description or rheology)
    TP = 1350.0 #mantle potential temperature [oC], at the LAB

    @property
    def ccapacity(self):
        '''Heat capacity [J/kg/K]: another name of the setting heat_capacity.'''
        return self.heat_capacity

    @ccapacity.setter
    def ccapacity(self, value):
        self.heat_capacity = value

    #Accordion regime [Myr]
    ti_convergence = 25.0
    ti_quiescence = 75.0
    dt_rifting2 = 60.0

    #Velocity boundary conditions
    fac_air = 10.0e3 #no velocity in the sticky air above thickness_sa - fac_air [m]
    h_v_reference = 'lithosphere' #or 'lithospheric_mantle': reference thickness of the region with constant velocity

    #Plots: profiles of temperature and YSE
    min_profile_width = 0.05 #fraction of Lx
    profile_colors = ['tab:orange', 'tab:purple', 'tab:green', 'tab:blue', 'r', 'k'] #the thickest lithosphere is the last one

    def __init__(self, svg_file, dx=None, dz=None, asthenosphere_label='asthenosphere', seed_label='weak_seed',
                 sampling=100.0, min_thickness=None, verbose=True, **settings):

        self.verbose = verbose
        self.svg_file = self._check_file('SVG file', svg_file)
        self._check_keys('setting', settings, EXPERIMENT_KEYS)
        setup = self._check_setup(read_setup(svg_file), svg_file, SETUP_LABEL)
        self.svg_version, self.inkscape_version = svg_versions(svg_file)
        self._check_svg_versions(svg_file, self.svg_version, self.inkscape_version, SVG_VERSION, INKSCAPE_VERSION)

        #domain and resolution
        self.Lx, self.Lz = setup['Lx'], setup['Lz']
        self.Nx = points_from_spacing(self.Lx, dx * 1.0e3, 'dx') if(dx) else setup['Nx']
        self.Nz = points_from_spacing(self.Lz, dz * 1.0e3, 'dz') if(dz) else setup['Nz']
        self._check_domain(self.Lx, self.Lz, self.Nx, self.Nz)
        self.x = np.linspace(0, self.Lx, self.Nx)
        self.z = np.linspace(self.Lz, 0, self.Nz)
        self.X, self.Z = np.meshgrid(self.x, self.z)
        self.dx = self.Lx / (self.Nx - 1)
        self.dz = self.Lz / (self.Nz - 1)

        #numerical experiment: argument > SVG > default
        self._experiment = setup.get('experiment', {})
        self.source = {}
        self._read_settings(settings)

        #layers
        self.asthenosphere_label = asthenosphere_label
        self.seed_label = seed_label
        self.min_thickness = min_thickness * 1.0e3 if(min_thickness is not None) else setup.get('min_thickness', MIN_THICKNESS)
        self.layers = layers_from_svg(svg_file, self.Lx, self.Lz, self.Nx, dz=sampling, setup_id=setup['id'],
                                      min_thickness=self.min_thickness, verbose=verbose)
        self._check_layers(self.layers, svg_file, asthenosphere_label)

        #water and magmatic layers: chosen by the user (label, or ID of the layer from the bottom), never guessed
        water = self._setting('water_layer', settings.get('water_layer'), None)
        if(water is not None and not str(water).strip().lstrip('-').isdigit() and not layer_indices(self.layers, str(water).strip())):
            #water not drawn in the domain (e.g. a template outside the setup shape): zero thickness below the sticky air
            self._insert_unit(str(water).strip(), len(self.layers) - 1)
        self.water_layer = self._layer_index('water_layer', water)
        if(self.water_layer is not None and self.water_layer != len(self.layers) - 2):
            self._fail(f"water_layer '{self.layers[self.water_layer].layer_label}' must be the unit right below the sticky air.")
        magmatic = self._setting('magmatic_layer', settings.get('magmatic_layer'), None)
        if(magmatic is not None and not str(magmatic).strip().lstrip('-').isdigit() and not layer_indices(self.layers, str(magmatic).strip())):
            #magmatic layer not drawn in the domain (a template outside the setup shape, e.g. a basalt with the
            #material of the SVG): zero thickness at the top of the solid layers, where the default basalt goes
            if(self._insert_unit(str(magmatic).strip(), len(self.layers) - 1 if(self.water_layer is None) else self.water_layer) is not None
               and self.water_layer is not None):
                self.water_layer += 1
        magmatic_layer = self._layer_index('magmatic_layer', magmatic)
        if(magmatic_layer is None and self.magmatism == 'on' and self.magmatism_extraction == 'on'):
            #Basalt layer with zero thickness at the top of the solid layers (below the water, if there is water,
            #otherwise below the sticky air): it can not be drawn in the SVG and receives the melt extracted by Mandyoc
            magmatic_layer = self._insert_basalt(len(self.layers) - 1 if(self.water_layer is None) else self.water_layer)
            if(self.water_layer is not None):
                self.water_layer += 1
        #magmatic_layer of param.txt: by default the top of the solid layers (below the water or the sticky air)
        top_solid = len(self.layers) - 2 if(self.water_layer is None) else self.water_layer - 1
        self.magmatic_layer = magmatic_layer if(magmatic_layer is not None) else top_solid

        #layers of sediments deposited by Mandyoc: zero thickness at the top of the solid layers (SP_vc.txt)
        self.sediment_layers = []
        if(self.sedimentation):
            self.sediment_layers = self._insert_sediments(top_solid + 1)
            if(self.sp_format == 'ggciag' and self.water_layer is not None):
                self._warn(f"water_layer '{self.layers[self.water_layer].layer_label}': the Mandyoc of ggciag has no water "
                           f"layer (water_interface is not written).")

        #sticky air and lithosphere along x
        #surface: top of the solid layers (base of the sticky air, or of the water); its interface ID is the
        #ID of the layer below it (interface i is the top of layer i)
        self.surface_id = len(self.layers) - 2 if(self.water_layer is None) else self.water_layer - 1
        self.thickness_sa = self.layers[self.surface_id + 1].interface #depth of the surface
        self.idx_ast = max(layer_indices(self.layers, asthenosphere_label))
        ## NOTE (Claude): the LAB was the top of the highest asthenosphere slot (idx_ast) in every column. In a
        ## deformed state (e.g. an SVG made from the final state of a simulation) the asthenosphere may get more
        ## than one slot, and where the highest one has zero thickness the LAB fell on the base of the crust (32 km
        ## instead of 129 km at the edges of the final state of the rift). The LAB is now found column by column
        ## (_find_lab); with one asthenosphere slot it is the same as before.
        self.idx_lab = self._find_lab()
        columns = np.arange(self.Nx)
        self.thickness_litho = np.zeros(self.Nx) #surface to the top of the asthenosphere (LAB)
        self.thickness_mantle = np.zeros(self.Nx) #lithospheric mantle (olivine rheologies above the asthenosphere)
        for i in range(len(self.layers) - 1):
            ## NOTE (Claude): top of the LAB slot of each column, instead of the top of idx_ast in all of them
            above = self.idx_lab == i
            self.thickness_litho[above] = self.layers[i+1].interface[above] - self.thickness_sa[above]
            ## NOTE (Claude): mantle above the LAB of each column, instead of above idx_ast in all of them
            if(self.layers[i].rheology_name in MANTLE_RHEOLOGIES):
                self.thickness_mantle += np.where(self.idx_lab < i, slot_thickness(self.layers, i, columns), 0.0)

        #position of the vertical profile of the scenario infos and the profiles of the plots
        self.x_profile = self.Lx / 4.0
        self.idx_profile = find_nearest(self.x, self.x_profile)
        self.profiles = self._find_profiles()

        #generic builder of tapioca (tapioca.inputs.MandyocBuilder) with the interfaces and the material fields of
        #the SVG; the temperature, the velocity and the YSE are built on it and the input files are exported by it
        self.builder = layers_to_builder(self.layers, self.Lx, self.Lz, Nz=self.Nz, seed_label=seed_label, verbose=False)
        self.builder.create_materials_fields()
        self.temperature_builder, self.velocity_builder = None, None

        self.T, self.VX, self.VZ = None, None, None

    ###########################################################################
    # Settings
    ###########################################################################

    def _setting(self, name, value, default):
        '''value (argument) > SVG > default'''
        if(value is not None):
            self.source[name] = 'argument'
            return value
        if(name in self._experiment):
            self.source[name] = 'SVG'
            return self._experiment[name]
        self.source[name] = 'default'
        return default

    def _layer_index(self, name, value):
        '''Index (ID from the bottom) of the layer given by its label (the top slot of the unit) or its ID; None if value is None.'''
        if(value is None):
            return None
        value = str(value).strip()
        if(value.lstrip('-').isdigit()):
            index = int(value)
            if(not 0 <= index < len(self.layers) - 1):
                self._fail(f"{name} = {index}: the layers below the sticky air have IDs 0 to {len(self.layers) - 2}.")
            return index
        indices = layer_indices(self.layers, value)
        if(not indices):
            self._fail(f"{name} = {value}: no unit labelled '{value}' in {self.svg_file}.")
        return max(indices)

    def _insert_basalt(self, index):
        '''Insert a basalt layer of zero thickness at index of the stack (below the layer at index); return its index.'''
        basalt = MandyocLayer('basalt', BASALT, Nx=self.Nx, title='Basalt', color='xkcd:dark grey blue',
                              rho=2900.0, C=1.0, H=9.0e-11)
        basalt.interface = np.copy(self.layers[index].interface) #base of basalt = base of the layer above it
        basalt.base_depth = basalt.interface.max()
        self.layers.insert(index, basalt)
        return index

    def _insert_unit(self, label, index):
        '''
        Insert a layer of zero thickness with the material of the unit labelled label (drawn in the domain or as a
        template outside it, see reader.unit_layer) at index of the stack (below the layer at index); return it,
        or None if the SVG has no unit with that label.
        '''
        layer = unit_layer(self.svg_file, label, self.Nx)
        if(layer is None):
            return None
        layer.interface = np.copy(self.layers[index].interface) #base = base of the layer above it
        layer.base_depth = layer.interface.max()
        self.layers.insert(index, layer)
        return layer

    def _insert_sediments(self, index):
        '''
        Insert one layer of zero thickness of sediment_layer for each time of sediment_times at index of the stack
        (the top of the solid layers, below the water or the sticky air); return their indices (bottom to top).
        '''
        for _ in self.sediment_times:
            if(self._insert_unit(self.sediment_layer, index) is None):
                self._fail(f"sediment_layer = {self.sediment_layer}: no unit labelled '{self.sediment_layer}' in "
                           f"{self.svg_file}. Draw a shape with its description (outside the setup shape, if the "
                           f"sediments are not in the initial scenario).")
        if(self.water_layer is not None):
            self.water_layer += len(self.sediment_times)
        return list(range(index, index + len(self.sediment_times)))

    def sedimentation_schedule(self):
        '''
        [(time [Myr], ID of the active layer of sediments, sedimentation rate)] at each time at which the active layer
        or the rate changes: the layers start at sediment_times and the rate changes at the times of
        sedimentation_rates. When tf_sed is before time_max, the last line (at tf_sed, rate 0) stops the
        sedimentation. Empty without sedimentation.
        '''
        if(not self.sedimentation):
            return []
        changes = [t for t, _ in self.sedimentation_rates if self.ti_sed < t < self.tf_sed]
        schedule = []
        for t in sorted(set(self.sediment_times) | set(changes)):
            k = max(k for k, ts in enumerate(self.sediment_times) if ts <= t + 1.0e-9) #active layer
            schedule.append((t, self.sediment_layers[k], self.rate_at(t)))
        if(self.tf_sed < self.time_max / 1.0e6 - 1.0e-9):
            schedule.append((self.tf_sed, self.sediment_layers[-1], 0.0))
        return schedule

    def time_conditions(self):
        '''
        Files of conditions that change in time (tapioca.inputs.GenericTimeConditions), set in the builder and written
        by write_all: scale_bcv.txt (accordion) and the files of the surface processes in the format of sp_format:
          - SP_vc: SP_vc.txt (each change of sedimentation_schedule: ID of the active layer, rate and Ld) and
            base_level.txt without header;
          - ggciag: sediment_layers.txt (the time and the ID of each layer; the header has the initial active layer),
            sedimentation_rate.txt (each change of the rate; only when sp_mode deposits sediments, Mandyoc starts with
            rate 0) and base_level.txt, each with the number of lines first.
        '''
        conditions = []
        if(self.variable_bcv):
            conditions.append(GenericTimeConditions('scale_bcv', self._bcv_changes()))
        if(self.sedimentation and self.sp_format == 'SP_vc'):
            conditions.append(GenericTimeConditions('SP_vc', [(t, (i, rate, self.Ld))
                                                              for t, i, rate in self.sedimentation_schedule()]))
        if(self.sedimentation and self.sp_format == 'ggciag'):
            conditions.append(GenericTimeConditions('sediment_layers', list(zip(self.sediment_times, self.sediment_layers))))
            if(self.sp_mode in SEDIMENTATION_MODES):
                rates = []
                for t, _, rate in self.sedimentation_schedule():
                    if(not rates or rates[-1][1] != rate):
                        rates.append((t, rate))
                conditions.append(GenericTimeConditions('sedimentation_rate', rates))
        if(self.base_level is not None):
            conditions.append(GenericTimeConditions('base_level', [(0.0, self.base_level), (self.time_max/1.0e6, self.base_level)],
                                                    header='count' if(self.sp_format == 'ggciag') else None))
        return conditions

    def _time_condition_text(self, param):
        '''Content of the file of the time condition param, or None if the scenario does not have it.'''
        return next((c.text() for c in self.time_conditions() if c.param == param), None)

    def sp_vc(self):
        '''Content of SP_vc.txt (sp_format = SP_vc; None without sedimentation), see time_conditions.'''
        return self._time_condition_text('SP_vc')

    def sediment_layers_file(self):
        '''Content of sediment_layers.txt (sp_format = ggciag; None without sedimentation), see time_conditions.'''
        return self._time_condition_text('sediment_layers')

    def sedimentation_rate_file(self):
        '''Content of sedimentation_rate.txt (sp_format = ggciag, sp_mode in SEDIMENTATION_MODES), see time_conditions.'''
        return self._time_condition_text('sedimentation_rate')

    def base_level_file(self):
        '''Content of base_level.txt (None without base_level), see time_conditions.'''
        return self._time_condition_text('base_level')

    def sp_files(self):
        '''Files of the surface processes ({file name: content}) written by write_all, in the format of sp_format.'''
        return {c.fname: c.text() for c in self.time_conditions() if c.param != 'scale_bcv'}

    def rate_at(self, t):
        '''Sedimentation rate at time t [Myr] (sedimentation_rates is constant in steps).'''
        return [rate for ts, rate in self.sedimentation_rates if ts <= t + 1.0e-9][-1]

    def _parse_rates(self, value):
        '''
        sedimentation_rate as [(time [Myr], rate)], sorted by time: a number (constant rate from ti_sed), a text
        "t1: rate1, t2: rate2, ..." (times in Myr or with units, e.g. "0 Myr: 25, 20 Myr: 50") or a list of (time, rate)
        or a dict {time: rate} (times in Myr).
        '''
        where = f"Setting of {self.svg_file}"
        if(isinstance(value, dict)):
            pairs = list(value.items())
        elif(isinstance(value, (list, tuple))):
            pairs = list(value)
        elif(isinstance(value, str) and ':' in value):
            pairs = []
            for item in value.split(','):
                if(':' not in item):
                    self._fail(f"sedimentation_rate: write each step as time: rate (found '{item.strip()}').")
                t, rate = item.split(':', 1)
                pairs.append((parse_quantity(t.strip(), 'time', where, 'sedimentation_rate', default_unit='Myr') * 1.0e-6, rate))
        else:
            pairs = [(self.ti_sed, value)]
        if(any(not isinstance(pair, (list, tuple)) or len(pair) != 2 for pair in pairs)):
            self._fail(f"sedimentation_rate: give a number or pairs of time and rate (found {value!r}).")
        rates = sorted((self._check_number('time of sedimentation_rate', t, 0),
                        self._check_number('sedimentation_rate', rate, 0)) for t, rate in pairs)
        if(not rates):
            self._fail("sedimentation_rate: no rate given.")
        if(len({t for t, _ in rates}) < len(rates)):
            self._fail("sedimentation_rate: two rates at the same time.")
        if(rates[0][0] > self.ti_sed + 1.0e-9):
            self._fail(f"sedimentation_rate: the first rate starts at {rates[0][0]:g} Myr, after ti_sed = {self.ti_sed:g} Myr; "
                       f"give the rate at ti_sed.")
        return rates

    def _read_settings(self, settings):
        where = f"Setting of {self.svg_file}"
        def setting(name, default, quantity=None, unit=None, scale=1.0):
            '''
            Setting (argument > SVG > default). A text with a unit of tapioca.UNITS[quantity] (without unit: unit) is
            converted to the base unit of the quantity and multiplied by scale (e.g. yr -> Myr).
            '''
            value = self._setting(name, settings.get(name), default)
            if(quantity is not None and isinstance(value, str)):
                value = parse_quantity(value, quantity, where, name, default_unit=unit) * scale
            return value

        self.regime = self._check_choice('regime', setting('regime', 'rifting'), REGIMES)
        self.velocity = self._check_number('velocity', setting('velocity', 1.0, 'velocity', 'cm/yr')) #cm/yr
        self.moving_boundary = self._check_choice('moving_boundary', setting('moving_boundary', 'right'),
                                                  ['left', 'right', 'symmetric'])
        time_max = setting('time_max', None, 'time', 'Myr', 1.0e-6) #Myr
        if(time_max is not None):
            time_max = self._check_number('time_max', time_max, 0, strict=True)
        self.quiescence_time = self._check_number('quiescence_time', setting('quiescence_time', 0.0, 'time', 'Myr', 1.0e-6), 0) #Myr

        self.magmatism = self._check_on_off('magmatism', setting('magmatism', 'off'))
        self.magmatism_extraction = self._check_on_off('magmatism_extraction', setting('magmatism_extraction', 'off'))
        self.magmatism_Ts0 = self._check_number('magmatism_Ts0', setting('magmatism_Ts0', 1080.0, 'temperature', 'oC')) #oC
        self._check_experiment(self.regime, self.moving_boundary, self.source, self.magmatism, self.magmatism_extraction)

        for name, minimum, strict in THERMAL_SETTINGS:
            setattr(self, name, self._check_number(name, setting(name, DEFAULT_PARAMS[name]), minimum, strict))

        preset = setting('preset_initial_temperature', None)
        if(preset is not None): #a relative path in the SVG is relative to the folder of the SVG
            base = os.path.dirname(os.path.abspath(self.svg_file)) if(self.source['preset_initial_temperature'] == 'SVG') else None
            preset = self._check_path('preset_initial_temperature', preset, base)
        self.preset_initial_temperature = preset

        custom = setting('custom_veloc_field', None)
        if(custom is not None): #a relative path in the SVG is relative to the folder of the SVG
            base = os.path.dirname(os.path.abspath(self.svg_file)) if(self.source['custom_veloc_field'] == 'SVG') else None
            custom = self._check_path('custom_veloc_field', custom, base)
            if(not os.path.isfile(custom) or not custom.endswith('.py')):
                self._fail(f"custom_veloc_field must be a Python script (.py): {custom}")
            if(self.regime == 'quiescence'):
                self._fail("custom_veloc_field: the quiescence regime has no velocity at the boundaries.")
        self.custom_veloc_field = custom

        self.temperature_solver = self._check_choice('temperature_solver', setting('temperature_solver', '2D'),
                                                     ['2D', 'columns'], aliases={'1D': 'columns'})

        #Velocity and time:
        #   rifting: extension, the moving boundary moves outwards with velocity
        #   compression: convergence, the moving boundary pushes the lithosphere inwards with velocity
        #   accordion (as in M300): the boundary velocity of the rifting regime is scaled in time by scale_bcv.txt
        #       0 - ti_convergence: extension
        #       ti_convergence - ti_quiescence: convergence (velocity x -1)
        #       ti_quiescence - tf_quiescence: quiescence (velocity x 0.01)
        #       tf_quiescence - time_max: extension again (dt_rifting2)
        #   quiescence: no velocity imposed at the boundaries (velocity_from_ascii = False)
        self.velocity_from_ascii = (self.regime != 'quiescence')
        self.variable_bcv = (self.regime == 'accordion')
        self.dt_quiescence = self.quiescence_time
        if(self.variable_bcv):
            self.tf_quiescence = self.ti_quiescence + self.dt_quiescence
            if(time_max is not None): #time_max sets the duration of the second extension
                self.dt_rifting2 = time_max - self.tf_quiescence
                if(self.dt_rifting2 <= 0):
                    self._fail(f"time_max = {time_max:g} Myr ends before the second extension starts ({self.tf_quiescence:g} Myr).")
            self.time_max = (self.tf_quiescence + self.dt_rifting2) * 1.0e6
        elif(time_max is not None):
            self.time_max = time_max * 1.0e6
        elif(self.regime == 'quiescence'):
            self.time_max = self.dt_quiescence * 1.0e6 if(self.dt_quiescence > 0) else 1000.0e6
        else: #rifting and compression
            self.time_max = 60.0e6

        self._read_surface_processes(settings, setting)

        self._print_verbose(f"Numerical experiment: {self.regime} ({self.source['regime']}), velocity = {self.velocity:g} cm/yr "
                    f"({self.source['velocity']}), moving boundary: {self.moving_boundary} ({self.source['moving_boundary']})")
        self._print_verbose(f"Magmatism: {self.magmatism} ({self.source['magmatism']}), extraction: {self.magmatism_extraction} "
                    f"({self.source['magmatism_extraction']}), Ts0 = {self.magmatism_Ts0:g} oC ({self.source['magmatism_Ts0']})")
        self._print_verbose("Thermal parameters: " + ', '.join(f"{name} = {getattr(self, name):g} ({self.source[name]})"
                                                               for name, _, _ in THERMAL_SETTINGS))
        if(self.preset_initial_temperature is not None):
            self._print_verbose(f"Initial temperature from the preliminary scenario {self.preset_initial_temperature} "
                        f"({self.source['preset_initial_temperature']})")
        if(self.custom_veloc_field is not None):
            self._print_verbose(f"Velocity field from {self.custom_veloc_field} ({self.source['custom_veloc_field']})")

    def _rates_text(self):
        '''sedimentation_rates as text: the rate, or "t1 Myr: rate1, t2 Myr: rate2, ..." when it changes.'''
        if(len(self.sedimentation_rates) == 1):
            return f"{self.sedimentation_rates[0][1]:g}"
        return ', '.join(f"{t:g} Myr: {rate:g}" for t, rate in self.sedimentation_rates)

    def _check_bool(self, name, value):
        '''True or False (also on/off, yes/no).'''
        if(isinstance(value, (bool, np.bool_))):
            return bool(value)
        key = str(value).strip().lower()
        if(key in ('true', 'on', 'yes')):
            return True
        if(key in ('false', 'off', 'no')):
            return False
        self._fail(f"{name} must be True or False (found '{value}').")

    def _read_surface_processes(self, settings, setting):
        '''
        Settings of the surface processes (after time_max): sp_surface_processes and sp_mode of param.txt, the
        layers of sediments (sediment_layer, ti_sed, tf_sed, dt_sed, sedimentation_rate and Ld, see sp_vc), the base
        level (base_level) and the parameters of SP_PARAMS, which are written to param.txt as given.
        '''
        self.sp_surface_processes = self._check_bool('sp_surface_processes', setting('sp_surface_processes', False))
        self.sp_mode = self._check_choice('sp_mode', setting('sp_mode', 'none'), SP_MODES)
        if(self.sp_surface_processes and self.sp_mode == 'none'):
            self._fail(f"sp_surface_processes = True needs sp_mode ({', '.join(SP_MODES[1:])}).")
        if(not self.sp_surface_processes and self.sp_mode != 'none'):
            self._warn(f"sp_mode = {self.sp_mode} is not used with sp_surface_processes = False.")
        #format of the files of the surface processes (default ggciag): SP_vc (SP_vc.txt, the Mandyoc of the sedimentation of
        #Theunissen, mandyoc_tests) or ggciag (sediment_layers.txt and sedimentation_rate.txt, the Mandyoc of ggciag)
        self.sp_format = self._check_choice('sp_format', setting('sp_format', 'ggciag'), SP_FORMATS)
        #layers of sediments: the modes that deposit sediments; with ggciag, every mode (Mandyoc reads
        #sediment_layers.txt whenever sp_surface_processes is True)
        self.sedimentation = self.sp_surface_processes and (self.sp_mode in SEDIMENTATION_MODES or self.sp_format == 'ggciag')

        sed_keys = ('sediment_layer', 'ti_sed', 'tf_sed', 'dt_sed', 'sedimentation_rate', 'Ld')
        self.sediment_layer = setting('sediment_layer', None)
        ti = self._check_number('ti_sed', setting('ti_sed', 0.0, 'time', 'Myr', 1.0e-6), 0) #Myr
        tf = setting('tf_sed', None, 'time', 'Myr', 1.0e-6)
        tf = self.time_max / 1.0e6 if(tf is None) else self._check_number('tf_sed', tf, ti, strict=True) #Myr
        dt = self._check_number('dt_sed', setting('dt_sed', 1.0, 'time', 'Myr', 1.0e-6), 0, strict=True) #Myr
        self.ti_sed, self.tf_sed, self.dt_sed = ti, tf, dt
        self.sedimentation_rates = self._parse_rates(setting('sedimentation_rate', 25.0)) #[(Myr, rate)]
        self.sedimentation_rate = self.rate_at(ti) #rate at ti_sed
        self.Ld = self._check_number('Ld', setting('Ld', 35.0e3, 'length', 'm'), 0, strict=True) #m
        base_level = setting('base_level', None, 'length', 'm')
        self.base_level = None if(base_level is None) else self._check_number('base_level', base_level) #m
        self.sp_params = {name: value for name in SP_PARAMS if (value := setting(name, None)) is not None}

        self.sediment_times = []
        if(self.sedimentation):
            if(self.sediment_layer is None):
                self._fail(f"sp_mode = {self.sp_mode} deposits layers of sediments: give sediment_layer, the label of "
                           f"the unit of the sediments (a shape that may be drawn outside the setup shape).")
            self.sediment_layer = str(self.sediment_layer).strip()
            if(tf > self.time_max / 1.0e6 + 1.0e-9):
                self._warn(f"tf_sed = {tf:g} Myr is after time_max = {self.time_max/1.0e6:g} Myr.")
            n = int(np.floor((tf - ti) / dt + 1.0e-9))
            if(n < 1):
                self._fail(f"No layer of sediments between ti_sed = {ti:g} and tf_sed = {tf:g} Myr with dt_sed = {dt:g} Myr.")
            self.sediment_times = [ti + k * dt for k in range(n)]
            if(self.sp_format == 'ggciag' and not np.isclose(self.Ld, 35.0e3)):
                self._warn(f"Ld = {self.Ld/1.0e3:g} km is not used: the Mandyoc of ggciag has Ld = 35 km in its code.")
            self._print_verbose(f"Sedimentation ({self.sp_mode}): {n} layers of {self.sediment_layer}, one every "
                                f"{dt:g} Myr from {ti:g} to {tf:g} Myr; sedimentation_rate = {self._rates_text()}, "
                                f"Ld = {self.Ld/1.0e3:g} km")
        else:
            unused = [key for key in sed_keys if self.source[key] != 'default']
            if(unused):
                self._warn(f"{', '.join(unused)} not used: the layers of sediments need sp_surface_processes = True "
                           f"and sp_mode in {', '.join(SEDIMENTATION_MODES)}.")

    ###########################################################################
    # Lithosphere-asthenosphere boundary
    ###########################################################################

    def _find_lab(self):
        '''
        Slot of the stack whose top is the LAB in each column (Nx): going down from the surface, the first
        asthenosphere slot with some thickness in the column. With one asthenosphere slot (a drawn initial
        scenario) it is that slot in every column. In a deformed state (e.g. the final state of a simulation)
        the asthenosphere may get more than one slot, with zero thickness in most columns, and the highest one
        would put the LAB at the base of the crust there; a lithospheric body inside the asthenosphere (e.g. a
        drip) stays below the LAB. Columns without asthenosphere take the highest asthenosphere slot.
        '''
        idx_lab = np.full(self.Nx, self.idx_ast)
        found = np.zeros(self.Nx, dtype=bool)
        for i in sorted(layer_indices(self.layers, self.asthenosphere_label), reverse=True):
            filled = ~found & (slot_thickness(self.layers, i, np.arange(self.Nx)) > 0)
            idx_lab[filled] = i
            found |= filled
        return idx_lab

    ###########################################################################
    # Profiles for the plots
    ###########################################################################

    ## I think its a very complex function for a simple task, all data is stored inside the fields dataset
    ## All these profiles could be obtained from the: MandyocBuilder.DTree.fields.ds
    def _find_profiles(self):
        '''
        Temperature and YSE profiles of the plots: one for each lithospheric domain (e.g. oceanic,
        non-cratonic, cratonic), i.e. each region with the same lithosphere thickness and the same unit at its
        base over at least min_profile_width (outside the weak seed). Each profile is at the centre of the widest
        region and is labelled with the title of the lithospheric mantle unit above the asthenosphere.
        Return a list of (column, label, color).
        '''
        layers, Nx = self.layers, self.Nx
        litho_km = np.round(self.thickness_litho / 1.0e3, 1)
        has_seed = np.zeros(Nx, dtype=bool)
        for i in layer_indices(layers, self.seed_label):
            has_seed |= slot_thickness(layers, i, slice(None)) > 0
        ## NOTE (Claude): unit at the base of the lithosphere above the LAB of each column (was above idx_ast)
        base_unit = np.array([next((layers[i].layer_label for i in range(self.idx_lab[ix] + 1, len(layers))
                                    if slot_thickness(layers, i, ix) > 0), '') for ix in range(Nx)])
        litho_key = np.array([f"{litho_km[ix]}|{base_unit[ix]}" if(not has_seed[ix]) else '' for ix in range(Nx)])
        edges = np.r_[0, np.nonzero(litho_key[1:] != litho_key[:-1])[0] + 1, Nx]
        regions = {} #(lithosphere thickness, unit at its base): (first column, last column + 1)
        for begin, end in zip(edges[:-1], edges[1:]):
            key = (litho_km[begin], base_unit[begin])
            if(litho_key[begin] != '' and (end - begin) * self.dx >= self.min_profile_width * self.Lx):
                if((key not in regions) or (end - begin > np.diff(regions[key])[0])):
                    regions[key] = (begin, end)

        profiles = []
        colors = self.profile_colors[-len(regions):] if(len(regions) <= len(self.profile_colors)) else self.profile_colors
        for k, ((thickness, unit), (begin, end)) in enumerate(sorted(regions.items())):
            ix = (begin + end - 1) // 2
            ## NOTE (Claude): mantle unit above the LAB of the column of the profile (was above idx_ast)
            i_mantle = next((i for i in range(self.idx_lab[ix] + 1, len(layers)) if slot_thickness(layers, i, ix) > 0), self.idx_lab[ix] + 1)
            name = layers[i_mantle].title.replace(' lith. mantle', '').replace(' lithospheric mantle', '')
            profiles.append((ix, f"{name} ({thickness:.0f} km)", colors[k % len(colors)]))
        if(len(profiles) == 0):
            profiles = [(self.idx_profile, f"x = {self.x_profile/1.0e3:.0f} km", 'r')]
        return profiles

    ###########################################################################
    # Initial temperature
    ###########################################################################

    def build_temperature(self):
        '''
        Initial temperature (Nz, Nx) built with a TemperatureFieldBuilder (self.temperature_builder), whose
        lithosphere goes from the surface (base of the sticky air) to the top of the asthenosphere (LAB):

        - without preset_initial_temperature: steady-state conduction in the lithosphere, with the
          conductivity and the radiogenic heat of each layer, 0 oC at the surface and TP at the LAB, and the
          mantle adiabat below the LAB (TemperatureFieldBuilder.apply_linear_temperature). temperature_solver
          '2D' (default) solves the whole lithosphere at once (solve_heat_2D), so lateral conduction follows the
          geometry of the layers; 'columns' solves each column alone (solve_heat_1D);
        - with preset_initial_temperature: the final thermal state of the preliminary scenario adapted to the
          lithosphere of each column.

        The field is also the temperature of the builder (self.builder.DTree.fields.temperature).
        '''
        ## NOTE (Claude): the LAB is given as its z along x (column by column, _find_lab) instead of the interface
        ## idx_ast, which is the LAB only when the asthenosphere has one slot
        self.temperature_builder = TemperatureFieldBuilder(self.builder, LAB_interface=-(self.thickness_sa + self.thickness_litho),
                                                           surface_interface=self.surface_id, c_cap=self.heat_capacity,
                                                           g=-self.gravity_acceleration, alpha=self.thermal_expansion_coefficient)
        if(self.preset_initial_temperature is None):
            self.temperature_builder.apply_linear_temperature(t_pot=self.TP)
            if(self.temperature_solver == '2D'):
                self.temperature_builder.solve_heat_2D()
            else:
                self.temperature_builder.solve_heat_1D()
        else:
            try:
                depth, T_preset, sa_preset, litho_preset, self.preset_source = read_preset_temperature(self.preset_initial_temperature)
            except FileNotFoundError as err:
                raise SVGInputError(f"preset_initial_temperature: {err}")
            #Without the interfaces.txt of the preliminary scenario, its lithosphere is the most common one of this scenario
            if(litho_preset is None):
                values, counts = np.unique(np.round(self.thickness_litho / 1.0e3), return_counts=True)
                litho_preset = values[np.argmax(counts)] * 1.0e3
                print(f"Warning: interfaces.txt of the preliminary scenario not found; its lithosphere is assumed to be "
                      f"{litho_preset/1.0e3:.0f} km thick (the most common in {self.svg_file}).")
            self.thickness_litho_preset = litho_preset
            self.temperature_builder.t_pot = self.TP
            #sticky air and lithosphere of each column from the surface and the LAB of the temperature builder
            tb = self.temperature_builder
            tb.temperature.values = adapt_preset_temperature(depth, T_preset, sa_preset, litho_preset, -tb.z,
                                                             -tb.surface_interface.values, tb.lithosphere_thickness.values,
                                                             self.thickness_mantle)
            self._print_verbose(f"Preset of initial temperature: {self.preset_source}")
            self._print_verbose(f"\tsticky air of the preliminary scenario: {sa_preset/1.0e3:.0f} km; reference lithosphere: "
                        f"{litho_preset/1.0e3:.0f} km (other lithosphere thicknesses are adapted)")
        self.builder.create_temperature_field(self.temperature_builder)
        self.T = self.builder.DTree.fields.temperature.values #(Nz, Nx)
        return self.T

    ###########################################################################
    # Velocity boundary conditions
    ###########################################################################

    def build_velocity(self):
        '''
        Initial velocity (Nz, Nx) with the boundary conditions of the regime (None in the quiescence regime),
        set by VelocityFieldBuilder.set_kinematic_regime (self.velocity_builder): constant velocity down to h_v_const
        below the surface (reference thickness + 20 km, where the lithosphere is thickest in rifting and accordion,
        as in M300, and at the moving boundary in compression, or at the thickest of both boundaries with
        moving_boundary = symmetric) and the return flow below it, so that the integral of the inflow is zero.
        Rifting and accordion are extensions; in compression the return flow leaves through the moving boundary
        and the opposite boundary is fixed. With moving_boundary = symmetric, both boundaries move with velocity/2
        (outwards in rifting and accordion, inwards in compression), each with its own return flow. In accordion,
        scale_bcv.txt changes the velocity in time.
        With custom_veloc_field, the script builds the velocity instead (see _build_custom_velocity).
        The velocity fields are those of the builder (self.builder.DTree.fields.vx and vz).
        '''
        if(not self.velocity_from_ascii):
            self.VX, self.VZ = None, None
            return None
        if(self.custom_veloc_field is not None):
            return self._build_custom_velocity()

        idx_bc = {'left': [0], 'right': [self.Nx - 1], 'symmetric': [0, self.Nx - 1]}[self.moving_boundary] #columns of the moving boundaries
        thickness_ref_x = self.thickness_litho if(self.h_v_reference == 'lithosphere') else self.thickness_mantle
        thickness_ref_bc = thickness_ref_x[idx_bc].max() if(self.regime == 'compression') else thickness_ref_x.max()
        self.h_v_const = thickness_ref_bc + 20.0e3 #thickness with constant velocity

        v = VelocityFieldBuilder(self.builder, units='cm/y')
        v.set_kinematic_regime('compression' if(self.regime == 'compression') else 'extension', self.velocity,
                              moving_boundary=self.moving_boundary, h_const=self.h_v_const, air_gap=self.fac_air,
                              surface=self.surface_id)
        self.builder.create_velocity_field(v)
        self.velocity_builder = v

        self.VX = self.builder.DTree.fields.vx.values
        self.VZ = self.builder.DTree.fields.vz.values
        return self.VX, self.VZ

    def _load_custom_velocity(self):
        '''Module of the script custom_veloc_field (loaded once).'''
        if(getattr(self, '_custom_velocity_module', None) is None):
            path = self.custom_veloc_field
            name = 'tapioca_custom_veloc_' + os.path.splitext(os.path.basename(path))[0]
            spec = importlib.util.spec_from_file_location(name, path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            if(not callable(getattr(module, 'build_velocity', None))):
                self._fail(f"custom_veloc_field: {path} has no function build_velocity(scenario).")
            self._custom_velocity_module = module
        return self._custom_velocity_module

    def _build_custom_velocity(self):
        '''
        Velocity field built by the script custom_veloc_field. The script defines

            build_velocity(scenario) -> VelocityFieldBuilder

        which sets the velocity at the boundaries of a VelocityFieldBuilder(scenario.builder, units='cm/y') (the
        net flux through the boundaries must be zero, e.g. with conservate_veloc) and returns it; the fields are
        then created by scenario.builder.create_velocity_field. Optionally, infos(scenario) returns a list of lines
        that describe the velocity field in infos_<folder>.txt.
        '''
        module = self._load_custom_velocity()
        v = module.build_velocity(self)
        if(not isinstance(v, VelocityFieldBuilder)):
            self._fail(f"custom_veloc_field: build_velocity of {self.custom_veloc_field} must return a VelocityFieldBuilder.")
        self.builder.create_velocity_field(v)
        self.velocity_builder = v
        self.h_v_const = None
        self.VX = self.builder.DTree.fields.vx.values
        self.VZ = self.builder.DTree.fields.vz.values
        net = v.integrate_normal_components()
        inflow = np.abs(self.VX[:, [0, -1]]).max() * self.Lz #scale of the flux through the sides
        if(abs(net) > 1.0e-6 * max(inflow, 1.0e-30)):
            print(f"Warning: the velocity field of {os.path.basename(self.custom_veloc_field)} does not conserve mass: "
                  f"net flux through the boundaries = {net:.3e} m2/s")
        return self.VX, self.VZ

    def scale_bcv(self):
        '''Content of scale_bcv.txt in the accordion regime (Mandyoc multiplies the boundary velocity by each factor at the given time).'''
        if(not self.variable_bcv):
            return None
        return self._time_condition_text('scale_bcv')

    def _bcv_changes(self):
        '''[(time [Myr], factor)] of scale_bcv.txt in the accordion regime.'''
        if(self.dt_quiescence > 0):
            return [(self.ti_convergence, -1.0), (self.ti_quiescence, 0.01), (self.tf_quiescence, -100.0),
                    (self.time_max/1.0e6, 0.01)]
        #no quiescence: convergence followed directly by extension
        return [(self.ti_convergence, -1.0), (self.ti_quiescence, -1.0), (self.time_max/1.0e6, 0.01)]

    def parameters(self, **params):
        '''
        Parameters of param.txt (dict): DEFAULT_PARAMS (tapioca.inputs.svg._variables) and those of the scenario
        (time_max, velocity, magmatism and thermal settings). The grid (nx, nz, lx, lz) and n_interfaces are set by
        MandyocBuilder.set_params. Any parameter of PARAMS_TYPES can be changed with a keyword argument,
        e.g. parameters(particles_per_element=100, rheology_model=9). This only changes the returned values: to
        change a thermal setting in the whole scenario, use write_all or the setting of the scenario.
        '''
        unknown = [key for key in params if key not in PARAMS_TYPES]
        if(unknown):
            raise SVGInputError(f"Unknown parameter(s) of param.txt: {', '.join(sorted(unknown))}")
        sp = {'sp_surface_processes': self.sp_surface_processes}
        if(self.sp_mode != 'none'):
            sp['sp_mode'] = self.sp_mode
        if(self.base_level is not None):
            sp['variable_base_level'] = True
        if(self.sedimentation and self.water_layer is not None and self.sp_format == 'SP_vc'):
            sp['water_interface'] = self.water_layer
        return {**DEFAULT_PARAMS,
                **sp,
                **self.sp_params,
                'time_max': self.time_max,
                'variable_bcv': self.variable_bcv,
                'velocity_from_ascii': self.velocity_from_ascii,
                'magmatism': self.magmatism,
                'magmatism_extraction': self.magmatism_extraction,
                'magmatism_Ts0': self.magmatism_Ts0,
                'magmatic_layer': self.magmatic_layer,
                **{name: getattr(self, name) for name, _, _ in THERMAL_SETTINGS},
                **params}

    def _update_thermal_settings(self, params):
        '''
        Thermal settings from the keyword arguments of write_all (params). They are also used by the temperature,
        so an initial temperature already built with other values can not be written with them.
        '''
        values = {name: self._check_number(name, params.get(name, getattr(self, name)), minimum, strict)
                  for name, minimum, strict in THERMAL_SETTINGS}
        builder = getattr(self, 'temperature_builder', None)
        if(self.T is not None and builder is not None):
            built = {'heat_capacity': builder.c_cap, 'thermal_expansion_coefficient': builder.alpha,
                     'gravity_acceleration': abs(builder.g)}
            for name, value in values.items():
                if(not np.isclose(value, built[name])):
                    self._fail(f"{name} = {value:g}, but the initial temperature was built with {built[name]:g}. "
                               f"Set {name} before build_temperature (or before the plots, which build it).")
        for name, value in values.items():
            setattr(self, name, value)
            if(name in params):
                self.source[name] = 'argument'

    ###########################################################################
    # Writing the files
    ###########################################################################

    def write_scripts(self, output_dir='.', system=None, ncores=None, mandyoc_options='',
                      petsc_dir=None, petsc_arch=None, mandyoc=None, post_processing=True,
                      tapioca_dir=None, post_processing_dir=None, petsc_env=PETSC_ENV, job_name=None, mail_user=None):
        '''
        Write the scripts of the scenario in output_dir:
          - run_linux.sh or run_mac.sh: runs Mandyoc with mpirun of PETSc (run_mac.sh first looks for the PETSc
            of the conda environment petsc_env), or run_hypatia.sh: a SLURM job on Hypatia, the cluster of the
            group (sbatch run_hypatia.sh), that runs Mandyoc with the mpirun of the modules of HYPATIA and
            activates the conda environment of HYPATIA for the post-processing. Then they move the outputs to one
            folder per variable with 0_organize_outputs.sh (on the Mac, caffeinate keeps it awake). With
            post_processing, it then creates the .nc files with src/postrunning/meshes2netcdf.jl and
            particles2netcdf.jl (last step as reference step), makes the snapshots and the video of the scenario
            with plot_snapshots.py and frames_generator.py of post_processing
            (https://github.com/jopmas/post_processing), zips the scenario with zipper.sh and removes the output
            files (only if the zip is fine and has all of them); a step is skipped if what it needs is not found;
          - 0_organize_outputs.sh: copy of src/postrunning/0_organize_outputs.sh;
          - zipper.sh: zips the scenario in <folder>.zip, skipping the files already in the zip.

        Parameters
        ----------
        system : str, optional
            'linux', 'mac' or 'hypatia'. Default: the system that runs tapioca ('linux' or 'mac').
        ncores : int, optional
            Number of MPI processes. Default: number of cores of this computer (hypatia: HYPATIA['ncores'],
            and SLURM_NTASKS when the job runs).
        mandyoc_options : str, optional
            Command line options of Mandyoc (e.g. '-seed 0,5,8 -strain_seed 0.0,1.0,1.0').
        petsc_dir, petsc_arch, mandyoc : str, optional
            PETSc folder and architecture (mpirun is $PETSC_DIR/$PETSC_ARCH/bin/mpirun) and the Mandyoc
            executable. Default: PETSC_DIR and PETSC_ARCH of the environment, or $HOME/opt/petsc and the
            only arch-* folder in it with mpirun (arch-label-optimized if there is none or more than one),
            and $HOME/opt/mandyoc/bin/mandyoc. Values set in the environment when the script runs override them.
        petsc_env : str or None, optional
            Conda environment whose PETSc run_mac.sh uses when PETSC_DIR is not set in the environment
            (default: 'petsc-env', the PETSc that Mandyoc is linked to; PETSC_ENV in the environment overrides it).
            The PETSc of petsc_dir and petsc_arch is used if the environment is not found. Ignored on Linux and
            when petsc_dir is given; None turns the lookup off.
        post_processing : bool, optional
            Add the post-processing to the run script (default: True).
        tapioca_dir, post_processing_dir : str, optional
            Folders of tapioca and of post_processing used by the post-processing. Default: TAPIOCA_DIR and
            POST_PROCESSING_DIR of the environment, or $HOME/opt/tapioca and $HOME/opt/post_processing. Values
            set in the environment when the script runs override them (also JULIA and PYTHON, the programs).
            On hypatia, the environment of this computer is not used: the defaults are in $HOME/opt of the
            cluster ($HOME/opt/petsc with HYPATIA['petsc_arch'], $HOME/opt/mandyoc/bin/mandyoc, ...).
        job_name, mail_user : str, optional
            hypatia only: name of the SLURM job (default: svg-<folder>) and e-mail that receives the messages of
            the job (default: TAPIOCA_MAIL_USER of the environment, or no messages).
        '''
        this_system = {'Linux': 'linux', 'Darwin': 'mac'}.get(platform.system())
        if(system is None):
            if(this_system is None):
                raise SVGInputError(f"Run scripts are written for {' or '.join(SYSTEMS)}; give system "
                                    f"(this computer is {platform.system()}).")
            system = this_system
        if(system not in SYSTEMS):
            raise SVGInputError(f"Unknown system '{system}'; use {' or '.join(SYSTEMS)}.")
        if(not os.path.isfile(ORGANIZE_SCRIPT)):
            raise SVGInputError(f"0_organize_outputs.sh not found in {os.path.dirname(os.path.abspath(ORGANIZE_SCRIPT))}.")

        if(system != 'mac' or petsc_dir is not None):
            petsc_env = None  #an explicit petsc_dir wins over the conda environment
        if(system == 'hypatia'):  #paths of the cluster, not of the environment of this computer
            petsc_dir = '$HOME/opt/petsc' if(petsc_dir is None) else petsc_dir
            petsc_arch = HYPATIA['petsc_arch'] if(petsc_arch is None) else petsc_arch
            tapioca_dir = '$HOME/opt/tapioca' if(tapioca_dir is None) else tapioca_dir
            post_processing_dir = '$HOME/opt/post_processing' if(post_processing_dir is None) else post_processing_dir
            ncores = HYPATIA['ncores'] if(ncores is None) else ncores
        if(petsc_dir is None):
            petsc_dir = os.environ.get('PETSC_DIR', '$HOME/opt/petsc')
        if(petsc_arch is None):
            petsc_arch = os.environ.get('PETSC_ARCH')
        if(petsc_arch is None):
            #the arch of $HOME/opt/petsc, only if the script is for this computer
            archs = glob.glob(os.path.join(os.path.expanduser('~'), 'opt', 'petsc', 'arch-*', 'bin', 'mpirun'))
            petsc_arch = (archs[0].split(os.sep)[-3] if(system == this_system and len(archs) == 1)
                          else 'arch-label-optimized')
        if(mandyoc is None):
            mandyoc = '$HOME/opt/mandyoc/bin/mandyoc'
        if(ncores is None):
            ncores = os.cpu_count() or 1
        if(tapioca_dir is None):
            tapioca_dir = os.environ.get('TAPIOCA_DIR', '$HOME/opt/tapioca')
        if(post_processing_dir is None):
            post_processing_dir = os.environ.get('POST_PROCESSING_DIR', '$HOME/opt/post_processing')

        os.makedirs(output_dir, exist_ok=True)
        path = lambda fname: os.path.join(output_dir, fname)

        run = HYPATIA_SCRIPT if(system == 'hypatia') else RUN_SCRIPT
        run = run.replace('@POST_PROCESSING@', POST_PROCESSING_SCRIPT if(post_processing) else NO_POST_PROCESSING_SCRIPT)
        run = run.replace('@PETSC_ENV_LOOKUP@', PETSC_ENV_LOOKUP if(petsc_env) else '')
        if(system == 'hypatia'):
            if(mail_user is None):
                mail_user = os.environ.get('TAPIOCA_MAIL_USER')
            if(job_name is None):
                job_name = f"svg-{os.path.basename(os.path.abspath(output_dir))}"
            sbatch_mail = (f"#SBATCH --mail-type=BEGIN,END,FAIL\n#SBATCH --mail-user={mail_user}\n" if(mail_user) else '')
            run = run.replace('@CONDA@', HYPATIA_CONDA_SCRIPT if(post_processing) else '')
            for key, value in [('SBATCH_MAIL', sbatch_mail), ('QOS', HYPATIA['qos']), ('TIME', HYPATIA['time']),
                               ('JOB_NAME', job_name), ('CONDA_ENV', HYPATIA['conda_env']),
                               ('MODULES', '\n'.join(f'module load {module}' for module in HYPATIA['modules']))]:
                run = run.replace(f'@{key}@', str(value))
        for key, value in [('PETSC_DIR', petsc_dir), ('PETSC_ARCH', petsc_arch), ('MANDYOC', mandyoc),
                           ('NCORES', int(ncores)), ('MANDYOC_OPTIONS', mandyoc_options),
                           ('LAUNCHER', RUN_LAUNCHER.get(system, '')), ('TAPIOCA_DIR', tapioca_dir),
                           ('POST_PROCESSING_DIR', post_processing_dir), ('PETSC_ENV', petsc_env)]:
            run = run.replace(f'@{key}@', str(value))
        scripts = {f'run_{system}.sh': run, 'zipper.sh': ZIP_SCRIPT}
        with open(ORGANIZE_SCRIPT) as f:
            scripts['0_organize_outputs.sh'] = f.read()

        for fname, text in scripts.items():
            with open(path(fname), 'w') as f:
                f.write(text)
            os.chmod(path(fname), 0o755)
        if(self.verbose):
            print(f"Scripts: {', '.join(scripts)} ({system}, {int(ncores)} cores"
                  f"{', with post-processing' if(post_processing) else ''})")

    def write_all(self, output_dir='.', plots=True, copy_svg=True, scripts=True, system=None, ncores=None,
                  mandyoc_options='', post_processing=True, zip_inputs=True, job_name=None, mail_user=None, **params):
        '''
        Write the input files of Mandyoc in output_dir: interfaces.txt, param.txt, input_temperature_0.txt,
        input_velocity_0.txt, vel_bc.txt, velz_bc.txt and scale_bcv.txt (accordion), SP_vc.txt (sedimentation) and
        base_level.txt (base_level), the scenario infos
        (infos_<folder>.txt), a copy of the SVG and the figures numerical_setup.png and initial_temperature.png.
        With scripts, also the run script of Mandyoc for system ('linux', 'mac' or 'hypatia'), with the post-processing
        if post_processing, 0_organize_outputs.sh and zipper.sh (see write_scripts).
        With zip_inputs, the files written are zipped in inputs_<folder>.zip (see zip_inputs).
        Keyword arguments change parameters of param.txt (see parameters); heat_capacity,
        thermal_expansion_coefficient and gravity_acceleration also change the initial temperature, so they
        must be given before it is built.
        '''
        self._update_thermal_settings(params)
        os.makedirs(output_dir, exist_ok=True)
        path = lambda fname: os.path.join(output_dir, fname)

        if(copy_svg):
            target = path(os.path.basename(self.svg_file))
            if(not os.path.exists(target) or not os.path.samefile(self.svg_file, target)):
                shutil.copy(self.svg_file, target)

        if(copy_svg and self.custom_veloc_field is not None): #the script of the velocity field goes with the SVG
            target = path(os.path.basename(self.custom_veloc_field))
            if(not os.path.exists(target) or not os.path.samefile(self.custom_veloc_field, target)):
                shutil.copy(self.custom_veloc_field, target)

        builder = self.builder
        builder.path = output_dir
        builder.export_interfaces(export='none')
        builder.set_params(**self.parameters(**params))
        builder.export_params()

        if(self.T is None):
            self.build_temperature()
        builder.export_field('temperature', header="T1\nT2\nT3\nT4")

        if(self.velocity_from_ascii):
            if(self.VX is None):
                self.build_velocity()
            #fields of the builder, which may have been changed by script after build_velocity
            self.VX, self.VZ = builder.DTree.fields.vx.values, builder.DTree.fields.vz.values
            X = self.X
            v0 = self.VX[(X == 0)]
            vf = self.VX[(X == self.Lx)]
            np.savetxt(path("vel_bc.txt"), np.array([v0, vf]).T)
            np.savetxt(path("velz_bc.txt"), self.VZ[(self.z == 0)].T)
            builder.export_field('velocity', header="v1\nv2\nv3\nv4")

        #files of conditions that change in time (scale_bcv.txt, files of the surface processes), written by the
        #builder with the others set in it by script (builder.set_time_conditions)
        for c in self.time_conditions():
            builder.set_time_conditions(c.param, c)
        builder.export_time_conditions()

        if(scripts):
            self.write_scripts(output_dir, system=system, ncores=ncores, mandyoc_options=mandyoc_options,
                               post_processing=post_processing, job_name=job_name, mail_user=mail_user)

        if(plots):
            self.plot_temperature(path('initial_temperature.png'))
            self.plot_setup(path('numerical_setup.png'))

        name = os.path.basename(os.path.abspath(output_dir))
        infos = self.infos(name)
        np.savetxt(path(f'infos_{name}.txt'), infos, fmt="%s")
        if(self.verbose):
            for line in infos:
                print(line)

        if(zip_inputs):
            self.zip_inputs(output_dir)

    #files written by write_all (those that exist are zipped by zip_inputs)
    INPUT_FILES = ['interfaces.txt', 'param.txt', 'input_temperature_0.txt', 'input_velocity_0.txt', 'vel_bc.txt',
                   'velz_bc.txt', 'scale_bcv.txt', 'SP_vc.txt', 'sediment_layers.txt', 'sedimentation_rate.txt',
                   'base_level.txt', 'numerical_setup.png', 'initial_temperature.png',
                   'run_linux.sh', 'run_mac.sh', 'run_hypatia.sh', '0_organize_outputs.sh', 'zipper.sh']

    def zip_inputs(self, output_dir='.'):
        '''
        Zip the input files of the scenario in output_dir/inputs_<folder>.zip: the files of write_all that are
        in output_dir (INPUT_FILES, infos_<folder>.txt and the copy of the SVG). It is apart from <folder>.zip,
        the zip of the inputs and outputs that zipper.sh makes after the run. Return the path of the zip.
        '''
        name = os.path.basename(os.path.abspath(output_dir))
        target = os.path.join(output_dir, f'inputs_{name}.zip')
        files = (self.INPUT_FILES + list(getattr(self.builder, 'time_conditions', {}))
                 + [f'infos_{name}.txt', os.path.basename(self.svg_file)])
        if(self.custom_veloc_field is not None):
            files.append(os.path.basename(self.custom_veloc_field))
        files = [fname for fname in dict.fromkeys(files) if os.path.isfile(os.path.join(output_dir, fname))]
        with zipfile.ZipFile(target, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
            for fname in files:
                zf.write(os.path.join(output_dir, fname), arcname=fname)
        if(self.verbose):
            print(f"Inputs zipped in {target} ({len(files)} files)")
        return target

    ###########################################################################
    # Scenario infos
    ###########################################################################

    ## This is a useful function that could be inside the MandyocBuilder class
    def infos(self, name=''):
        '''Description of the scenario (list of lines).'''
        layers, ix = self.layers, self.idx_profile
        lines = ['SCENARIO INFOS:', ' ', f'Name: {name}',
                 f"Built with tapioca {__version__} from {os.path.basename(self.svg_file)}: SVG {self.svg_version or '(not declared)'}, "
                 f"saved by Inkscape {self.inkscape_version.split()[0] if self.inkscape_version else '(not declared)'} "
                 f"(tested: SVG {SVG_VERSION}, Inkscape {INKSCAPE_VERSION})",
                 ' ', 'Domain parameters:',
                 f"\tLx: {self.Lx*1.0e-3} km", f"\tLz: {self.Lz*1.0e-3} km", f"\tNx: {self.Nx}", f"\tNz: {self.Nz}",
                 f"Resolution dx x dz: {1.0e-3*self.dx} x {1.0e-3*self.dz} km2", ' ',
                 f'Time limit: {self.time_max/1.0e6} Myr']
        if(self.variable_bcv):
            lines += [f"Extension: 0 - {self.ti_convergence} Myr",
                      f"Convergence: {self.ti_convergence} - {self.ti_quiescence} Myr",
                      f"Time of quiescence: {self.dt_quiescence} Myr ({self.ti_quiescence} - {self.tf_quiescence} Myr)",
                      f"Second extension: {self.dt_rifting2} Myr ({self.tf_quiescence} - {self.time_max/1.0e6} Myr)", ' ']
        lines.append(f'Layers (bottom to top) from {self.svg_file}, thickness at x = {self.x_profile/1.0e3:.0f} km:')
        sediments = set(self.sediment_layers)
        for i, layer in enumerate(layers):
            if(i in sediments):
                if(i == self.sediment_layers[0]):
                    lines.append(f"\t{i}-{self.sediment_layers[-1]}: {layer.layer_label} ({layer.rheology_name}): "
                                 f"{len(sediments)} layers of zero thickness deposited by Mandyoc "
                                 f"({'SP_vc.txt' if(self.sp_format == 'SP_vc') else 'sediment_layers.txt'}), "
                                 f"C = {layer.params['C']}, rho = {layer.params['rho']} kg/m3, H = {layer.params['H']:.4g} W/kg")
                continue
            lines.append(f"\t{i}: {layer.layer_label} ({layer.rheology_name}): h = {slot_thickness(layers, i, ix)/1.0e3:.1f} km, "
                         f"C = {layer.params['C']}, rho = {layer.params['rho']} kg/m3, H = {layer.params['H']:.4g} W/kg")
        lines.append(f"\tlithosphere: {self.thickness_litho[ix]/1.0e3:.1f} km")
        lines.append(f"\tmagmatic layer: {self.magmatic_layer} ({layers[self.magmatic_layer].layer_label})")
        lines.append("Temperature and YSE profiles in numerical_setup.png:")
        for ixp, label, color in self.profiles:
            lines.append(f"\t{label}: x = {self.x[ixp]/1.0e3:.0f} km")
        lines += [' ', f"Numerical experiment: {self.regime}", f"Velocity field: {self.velocity_from_ascii}"]
        if(self.velocity_from_ascii and self.VX is not None and self.custom_veloc_field is not None):
            lines.append(f"initial velocity: {self.velocity:g} cm/yr, moving boundary: {self.moving_boundary}")
            lines.append(f"\tcustom velocity field: {self.custom_veloc_field}")
            module = self._load_custom_velocity()
            if(callable(getattr(module, 'infos', None))):
                lines += [f"\t{line}" for line in module.infos(self)]
        elif(self.velocity_from_ascii and self.VX is not None):
            lines.append(f"initial velocity: {self.velocity:g} cm/yr")
            lines.append(f"\tmoving boundary: {self.moving_boundary}, constant velocity down to {self.h_v_const/1.0e3:.0f} km below the surface")
        lines += [f"Variable velocity field: {self.variable_bcv}", 'Initial temperature field setup:',
                  f"\tPreset of initial temperature field: {self.preset_initial_temperature is not None}"]
        if(self.preset_initial_temperature is None):
            lines += [f"\tSteady-state conduction in the lithosphere ({self.temperature_solver} solver)",
                      f"\tMantle potential temperature (at the LAB): {self.TP} oC"]
        elif(hasattr(self, 'preset_source')):
            lines += [f"\tPreliminary scenario: {self.preset_source}",
                      f"\tReference lithosphere of the preliminary scenario: {self.thickness_litho_preset/1.0e3:.0f} km"]
        lines += [f"\tHeat capacity: {self.heat_capacity:g} J/kg/K",
                  f"\tThermal expansion coefficient: {self.thermal_expansion_coefficient:g} 1/K",
                  f"\tGravity: {self.gravity_acceleration:g} m/s2"]
        lines += [' ', f'magmatism: {self.magmatism}, extraction: {self.magmatism_extraction}, Ts0: {self.magmatism_Ts0:g} oC']
        lines.append(f"Surface processes: {self.sp_surface_processes}, sp_mode: {self.sp_mode}"
                     + (f", files: {', '.join(self.sp_files())} (sp_format = {self.sp_format})" if(self.sp_files()) else ""))
        if(self.sedimentation):
            lines += [f"\tsediments: {len(self.sediment_layers)} layers of {self.sediment_layer} (IDs {self.sediment_layers[0]}-"
                      f"{self.sediment_layers[-1]}), one every {self.dt_sed:g} Myr from {self.ti_sed:g} to {self.tf_sed:g} Myr",
                      f"\tsedimentation_rate: {self._rates_text()}, Ld: {self.Ld/1.0e3:g} km"]
        if(self.water_layer is not None):
            lines.append(f"\twater layer: {self.water_layer} ({layers[self.water_layer].layer_label})")
        if(self.base_level is not None):
            lines.append(f"\tbase level: {self.base_level:g} m (base_level.txt)")
        if(self.sp_params):
            lines.append("\tparam.txt: " + ', '.join(f"{k} = {v}" for k, v in self.sp_params.items()))
        lines.append('Settings (argument, SVG or default): ' + ', '.join(f"{k}: {v}" for k, v in self.source.items()))
        return lines

    ###########################################################################
    # Plots
    ###########################################################################

    def yield_strength_envelope(self, ix, sr=1.0e-15):
        '''YSE [Pa] along the vertical profile at column ix (strain rate sr [1/s]), from MandyocBuilder.calculate_YSE.'''
        if(self.T is None):
            self.build_temperature()
        fields = self.builder.DTree.fields
        if('YSE' not in fields or fields.YSE.attrs.get('strain_rate') != sr):
            self.builder.calculate_YSE(sr=sr, g=-self.gravity_acceleration)
        return self.builder.DTree.fields.YSE.isel(x=ix).values

    ## It's a good function to be inside the TemperatureFieldBuilder
    ## This function should get **kwargs to set the plot
    ## TemperatureFieldBuilder has all interfaces encapsulated
    def plot_temperature(self, fname='initial_temperature.png'):
        '''Initial temperature field with the 600, 1000 and 1300 oC isotherms and the layer interfaces.'''
        T = self.T if(self.T is not None) else self.build_temperature()
        X, Z, x = self.X, self.Z, self.x
        sa = self.thickness_sa[self.idx_profile]

        fig, ax0 = plt.subplots(nrows=1, ncols=1, constrained_layout=True, figsize=(16, 6))
        im = ax0.contourf(X / 1.0e3, (-Z + sa) / 1.0e3, T, levels=np.arange(0, np.max(T) + 100, 100), cmap='inferno')
        isotherms = ax0.contour(X / 1.0e3, (-Z + sa) / 1.0e3, T, levels=[600, 1000, 1300], colors='w', linewidths=1.0, linestyles='--')
        ax0.clabel(isotherms, fmt='%d °C', fontsize=10)
        for layer in self.layers[1:]:
            ax0.plot(x / 1.0e3, (-layer.interface + sa) / 1.0e3, color='k', lw=0.6)
        ax0.set_xlim(0, self.Lx / 1.0e3)
        ax0.set_ylim(-self.Lz/1000 + sa/1000, 0 + sa/1000)
        ax0.tick_params(labelsize=12)
        ax0.set_xlabel("Distance [km]", fontsize=14)
        ax0.set_ylabel("Depth [km]", fontsize=14)
        ax0.set_title(f"Initial temperature field ({os.path.basename(self.svg_file)})", fontsize=14)
        cbar = fig.colorbar(im, orientation='vertical', ax=ax0, pad=0.01)
        cbar.set_label("Temperature [°C]", fontsize=12)
        cbar.ax.tick_params(labelsize=10)
        fig.savefig(fname, dpi=300)
        plt.close(fig)

    def plot_setup(self, fname='numerical_setup.png'):
        '''
        Numerical setup: layers with the colours of the SVG, legend with the properties, weak seed,
        temperature and YSE profiles of each lithospheric domain, friction angle and velocity boundary conditions.
        '''
        T = self.T if(self.T is not None) else self.build_temperature()
        if(self.velocity_from_ascii and self.VX is None):
            self.build_velocity()
        with plt.rc_context({'xtick.labelsize': 18, 'ytick.labelsize': 18}):
            self._plot_setup(fname, T)

    def _plot_setup(self, fname, T):
        layers, x, z, Lx, Lz, Nx = self.layers, self.x, self.z, self.Lx, self.Lz, self.Nx
        sa = self.thickness_sa[self.idx_profile]
        dx, dz = self.dx, self.dz

        fig, axs = plt.subplots(1, 1, figsize=(14, 6))
        ylimplot = [-Lz/1000 + sa/1000, 0 + sa/1000]

        for i, layer in enumerate(layers):
            top = layers[i+1].interface if i < len(layers)-1 else np.zeros(Nx)
            axs.fill_between(x/1000, (-layer.interface + sa)/1000, (-top + sa)/1000, color=layer.color, lw=0)

        axs.set_xlim(0, Lx/1000)
        axs.set_xticks([])
        axs.set_yticks([])
        axs.set_ylim(ylimplot)
        axs.set_xlabel(f"Distance = {Lx/1000:.0f} km\nResolution = {dx/1000:.0f} km", fontsize=16)
        axs.set_ylabel(f"Depth = {Lz/1000:.0f} km\nResolution = {dz/1000:.0f} km", fontsize=16)

        #legend: one entry per label, from top to bottom; thickness h:
        #   bodies (units that do not reach both sides of the domain, e.g. a basin): maximum thickness (depocentre)
        #   bottom unit (asthenosphere): range of thickness along x
        #   other layers: most common thickness along x (their reference thickness)
        legend_layers = {}
        for layer in layers[::-1]:
            if(layer.layer_label != self.seed_label):
                legend_layers.setdefault(layer.layer_label, layer)

        def legend_thickness(label):
            th = layer_thickness(layers, label, slice(None)) / 1.0E3
            if(th[0] <= 0 or th[-1] <= 0):
                return f"{th.max():.0f} km"
            if(label == layers[0].layer_label and np.round(th.min()) != np.round(th.max())):
                return f"{th.min():.0f}-{th.max():.0f} km"
            values, counts = np.unique(np.round(th[th > 0], 1), return_counts=True)
            return f"{values[np.argmax(counts)]:.0f} km" if(len(values)) else "0 km"

        legend_elements = []
        for label, layer in legend_layers.items():
            legend_label = (f"{layer.title}\n{layer.params['C']:g} x {layer.rheology_name.replace('_', ' ')}\n"
                            + fr"$\rho$ = {layer.params['rho']:.0f} kg/m³" + f"\n$h$ = {legend_thickness(label)}")
            legend_elements.append(Patch(facecolor=layer.color, edgecolor='black', label=legend_label))
        fig.subplots_adjust(bottom=0.25) #the legend is placed at the end (_legend_below), below the label of the domain


        ## Idea: get the weak seed from the MandyocBuilder.DTree.interfaces, where weakening_seed > 0
        #weak seed (one annotation per body, i.e. per group of contiguous columns)
        idx_seed = layer_indices(layers, self.seed_label)
        if(len(idx_seed) > 0):
            seed_bases = np.array([layers[i].interface for i in idx_seed])
            seed_tops = np.array([layers[i+1].interface if i < len(layers)-1 else np.zeros(Nx) for i in idx_seed])
            seed_present = seed_bases > seed_tops
            seed_top = np.where(seed_present, seed_tops, np.inf).min(axis=0)
            seed_base = np.where(seed_present, seed_bases, -np.inf).max(axis=0)
            seed_cols = seed_present.any(axis=0)
            seed_title = layers[idx_seed[0]].title.lower()
            edges = np.diff(np.r_[0, seed_cols.astype(int), 0])
            for begin, end in zip(np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0]):
                cols = slice(begin, end)
                x_seed = np.mean(x[cols])
                width_seed = (end - begin) * dx
                thickness_seed = np.max(seed_base[cols] - seed_top[cols])
                axs.text(x_seed/1000, -(np.max(seed_base[cols]) - sa)/1000 - 2,
                         f'{seed_title}\n{width_seed/1000:.0f}x{thickness_seed/1000:.0f} km²',
                         color='k', fontsize=12, ha='center', va='top')

        #temperature profiles
        axt = axs.inset_axes((0.205, 0, 0.18, 1))
        for ix, label, color in self.profiles:
            axt.plot(T[:, ix], (-z + sa) / 1.0e3, "-", color=color, label=label)
            axs.plot(x[ix] / 1.0e3, sa / 1.0e3 * 0.15, 'v', color=color, ms=6) #position of the profile
        axt.legend(loc='lower left', fontsize=8, framealpha=0.9)
        axt.grid(visible=True, axis='x', which='both', ls='--', color='red', alpha=0.3)
        axt.set_ylim(ylimplot)
        axt.set_yticks([])
        axt.set_xticks(np.linspace(0, 1800, 7))
        axt.patch.set_alpha(0)
        axt.xaxis.set_ticks_position('top')
        axt.tick_params(labelsize=8)
        axt.tick_params(axis='x', colors='k')
        axt.spines['left'].set_visible(False)
        axt.spines['right'].set_visible(False)
        axt.set_title('Temperature [°C]', color='k')

        #yield strength envelopes
        axsg = axs.inset_axes((0.605, 0, 0.13, 1))
        for ix, label, color in self.profiles:
            axsg.plot(self.yield_strength_envelope(ix)/1e9, -(z - sa)/1e3, '-', color=color, label=label)
        axsg.grid(visible=True, axis='x', which='both', ls='--', color='gray', alpha=0.8)
        axsg.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        axsg.set_xlim(-0.1, 1.1)
        axsg.set_ylim(ylimplot)
        axsg.set_yticks([])
        axsg.patch.set_alpha(0)
        axsg.xaxis.set_ticks_position('top')
        axsg.spines['left'].set_visible(False)
        axsg.spines['right'].set_visible(False)
        axsg.set_title(r'$\sigma_{YSE}$ [GPa]')
        axsg.tick_params(labelsize=8)

        #effective friction angle
        axsf = axs.inset_axes((0.800, 0.15, 0.10, 0.30)) #above the base, so its labels do not cover the velocity there
        axsf.patch.set_alpha(0.4)
        fisize = 10
        axsf.plot(np.array([0, 0.25, 0.75, 1.0]), np.array([15, 15, 2, 2]), 'k-')
        axsf.set_xlim([0.10, 0.9])
        axsf.set_ylim([0, 17])
        axsf.set_xticks([0.25, 0.5, 0.75])
        axsf.set_xticklabels([0.05, ' ', 1.05])
        axsf.tick_params('x', top=True, labeltop=False)
        axsf.xaxis.set_ticks_position('bottom')
        axsf.set_xlabel(r"$\varepsilon$", fontsize=fisize)
        axsf.set_yticks([2, 15])
        axsf.set_yticklabels(['2°', '15°'])
        axsf.set_ylabel(r'$\Phi_{\mathrm{eff}}$', fontsize=fisize)
        axsf.tick_params(labelsize=fisize)
        axsfC = axsf.twinx()
        axsfC.patch.set_visible(False) #the background of the twin axes would hide the curve
        axsfC.set_ylim([0, 17])
        axsfC.set_yticks([2, 15])
        axsfC.set_yticklabels([4, 20])
        axsfC.set_ylabel('Cohesion [MPa]', fontsize=fisize, rotation=90)
        axsfC.tick_params(labelsize=fisize)

        #velocity boundary conditions: insets centred on the edges, so v = 0 is on the boundary
        if(self.velocity_from_ascii):
            scale = 100.0 * SEC_PER_YEAR #m/s to cm/yr
            fac = 0.96
            for column, origin, title in ((-1, fac, f'v = {self.velocity} [cm/y]'), (0, -(1-fac), None)):
                vmax = np.round(max(abs(self.VX[:, column] * scale)), 0)
                vmax = vmax if(vmax > 0) else 1.0
                axv = axs.inset_axes((origin, 0, (1-fac)*2, 1))
                #only the hatch is drawn: the hatch takes the edge color (with alpha=0 it is also invisible in matplotlib >= 3.11)
                axv.fill_betweenx((-z + sa) / 1.0e3, self.VX[:, column] * scale, 0, facecolor='none', edgecolor='k', linewidth=0, hatch='---')
                axv.set_ylim(ylimplot)
                axv.set_yticks([])
                axv.set_xticks([])
                axv.patch.set_alpha(0)
                if(title):
                    axv.set_title(title, fontsize=8)
                axv.set_xlim([-vmax, vmax])
                for side in ['left', 'right', 'bottom'] + (['top'] if(title is None) else []):
                    axv.spines[side].set_visible(False)
            #the hatches of the left boundary go beyond the domain where VX < 0 (clipped by the inset): the y label
            #(its right side is at the given x) is 0.5% of the width to the left of them
            vmin_left = max(-self.VX[:, 0].min() * scale, 0.0)
            hatch_left = (1 - fac) * min(vmin_left / vmax, 1.0)
            axs.yaxis.set_label_coords(-hatch_left - 0.005, 0.5)

            #velocity across the base (e.g. outflow of a custom velocity field): vertical lines in an inset centred on the base
            vz_base = self.VZ[np.argmax(self.z), :] * scale #the base is the deepest row
            vz_max = np.abs(vz_base).max()
            if(vz_max > 0):
                axb = axs.inset_axes((0, -(1-fac), 1, (1-fac)*2))
                axb.set_zorder(axs.get_zorder() - 0.5) #below the other insets
                axb.fill_between(x / 1.0e3, vz_base, 0, facecolor='none', edgecolor='k', linewidth=0, hatch='|||')
                axb.set_xlim(0, Lx / 1.0e3)
                axb.set_ylim(-vz_max, vz_max)
                axb.set_xticks([])
                axb.set_yticks([])
                axb.patch.set_alpha(0)
                for side in axb.spines:
                    axb.spines[side].set_visible(False)
                vz_mean = vz_base[vz_base != 0].mean()
                axb.text(1.005, 0.25, f'vz = {vz_mean:.3g}\n[cm/y]', transform=axb.transAxes, fontsize=8, ha='left', va='center')
                axs.xaxis.labelpad = 15 #the label of the domain below the lines

        self._legend_below(fig, axs, legend_elements)
        fig.savefig(fname, bbox_inches="tight", dpi=300)
        plt.close(fig)

    @staticmethod
    def _legend_below(fig, axs, handles, pad=0.02):
        '''
        Legend of the layers below the label of the domain, in as few rows as fit in the width of the domain (the
        entries are spread evenly over the rows; one entry per row if even a single one is wider).
        '''
        fig.canvas.draw() #positions of the labels after the layout
        renderer = fig.canvas.get_renderer()
        width = axs.get_window_extent(renderer).width
        label_bottom = axs.xaxis.label.get_window_extent(renderer).y0 #pixels
        top = axs.transAxes.inverted().transform((0, label_bottom))[1] - pad #axes fraction
        for rows in range(1, len(handles) + 1):
            ncol = int(np.ceil(len(handles) / rows))
            legend = axs.legend(handles=handles, ncol=ncol, loc='upper center', bbox_to_anchor=(0.5, top),
                                frameon=False, title='Layers properties', title_fontsize=15, fontsize=8,
                                columnspacing=1.0)
            if(legend.get_window_extent(renderer).width <= width or ncol == 1):
                return legend
