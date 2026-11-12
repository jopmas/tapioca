# Scenarios from SVG files

One example of each numerical experiment of `tapioca.inputs.svg.SVGScenario`. Open the SVGs in Inkscape
(*Object > Object Properties*) to see the label, title and description of each unit, and the description of
the `setup` rectangle (in the `domain` layer) with the domain, the resolution and the experiment. The units
are in the `lithologies` layer. Every unit gives its reference thickness (`h = ... km`), except the
asthenosphere, which fills the rest of the domain, the layers whose thickness changes smoothly along x (e.g. the
thickened crust of `rifting_sedimentation`) and the templates drawn outside the domain.

| Folder | Scenario | Experiment |
|---|---|---|
| `rifting` | continental lithosphere (155 km) with a weak seed | rifting, 1 cm/yr, magmatism on |
| `rifting_symmetric` | the same lithosphere, with the weak seed at the centre | rifting, 1 cm/yr, both boundaries moving outwards at 0.5 cm/yr (`moving_boundary = symmetric`) |
| `rifting_with_basin` | continental lithosphere thinned below a sedimentary basin, with a weak seed (lithologies drawn at half horizontal scale) | rifting, 1 cm/yr, magmatism on |
| `compression` | oceanic lithosphere - continental margin ramp - non-cratonic lithosphere - craton, with a dipping weak zone | compression, the oceanic plate pushed from the left at 1 cm/yr |
| `accordion` | non-cratonic lithosphere between two cratons, with a weak seed (M300) | accordion (extension, convergence, quiescence and extension) |
| `quiescence` | cratonic keel with a mobile belt in a 4000 km domain | quiescence, 1000 Myr |
| `rifting_sedimentation` | continental lithosphere (120 km) with the crust 5 km thicker at the centre, as the scenario 0.25cm_35km_lc0.5 | symmetric rifting, 0.5 cm/yr, 70 Myr, sedimentation (`theunissen_sedimentation`) with 70 layers of sediments |

To create the input files of Mandyoc:

```bash
cd rifting
tapioca-svg rifting.svg -o inputs
```

The settings of the SVG can be changed on the command line, e.g. `tapioca-svg accordion.svg --quiescence-time 0 --dx 1 --dz 1`.

`moving_boundary = symmetric` (in the SVG, as in `rifting_symmetric`, or `--moving-boundary symmetric`)
puts the reference frame at the centre of the domain: both boundaries move with half of the velocity, outwards in
rifting and accordion and inwards in compression, e.g. `tapioca-svg compression.svg -o inputs --moving-boundary symmetric`.
In compression with `left` or `right`, the opposite boundary of the moving one is fixed.

The `rifting` and `quiescence` SVGs have the line `preset_initial_temperature = ../preliminary` commented
out: uncomment it (or use `--preset-initial-temperature`) with the folder of a preliminary scenario to start
from its final thermal state instead of the steady-state conduction in the lithosphere.

Sedimentation (`rifting_sedimentation`): with `sp_surface_processes = True` and `sp_mode = theunissen_sedimentation`
(or `sedimentation_rate_limited`), tapioca adds one layer of zero thickness of `sediment_layer` at the top of the solid
layers for each interval of `dt_sed` between `ti_sed` (default 0) and `tf_sed` (default `time_max`).
`sedimentation_rate` (m²/yr) may change in time, in steps: `sedimentation_rate = 0 Myr: 25, 20 Myr: 50` (each rate
from its time on). The sediments, the water and the extrusive basalt are drawn once, as templates below the domain (in
the `templates` layer): their shapes only give the material of the layers, so the 70 layers of sediments are not drawn.

The files depend on the version of Mandyoc (`sp_format`, or `--sp-format`):

- `ggciag` (default), the Mandyoc of ggciag (e.g. the one compiled in `~/opt/mandyoc`): `sediment_layers.txt`,
  `sedimentation_rate.txt` and `base_level.txt`, each with the number of lines first. Its `Ld` is fixed at 35 km and
  it has no water layer, so `water_interface` is not written.
- `SP_vc`, the Mandyoc of the group with the sedimentation of Theunissen (`mandyoc_tests`, on Hypatia): `SP_vc.txt`
  (time, ID of the active layer, `sedimentation_rate` and `Ld` per line), `base_level.txt` without the number of lines
  and `water_interface` in `param.txt`. `tapioca-svg rifting_sedimentation.svg --sp-format SP_vc --water-layer water`
  gives the inputs of the scenario 0.25cm_35km_lc0.5.

The other
parameters of the surface processes (`sp_d_c`, `strain_sed`, `continental_slope`, `a2l`, ...) may also be written in
the description of the `setup` rectangle; they go to `param.txt` as given.
