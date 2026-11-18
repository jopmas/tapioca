# Classical examples of Mandyoc as SVG scenarios

The examples of Mandyoc (`mandyoc/examples`, one `create_input*.ipynb` each) drawn as SVG scenarios of
`tapioca.inputs.svg.SVGScenario`. Each SVG has the geometry, the resolution (`Nx`, `Nz`) and the materials of its
notebook. The description of the `setup` shape also holds the parameters of `param.txt` of the notebook that differ
from the defaults of tapioca (e.g. `rheology_model = 9`). The comments at the top of each SVG describe the example.

The SVGs are written by `make_classical_examples.py`:

```bash
python make_classical_examples.py
```

| Folder | Example | Settings of tapioca it uses |
|---|---|---|
| `continental_rift` | crust 20 + 20 km, lithosphere 130 km, weak seed in the lithospheric mantle; rifting at 1 cm/yr, 100 Myr | |
| `continental_rift_with_surface_processes` | the same at 2 km, with surface processes | `sp_mode = diffusion` (the notebook uses the old `sp_mode = 5`) |
| `continental_rift_non_dimensional/model1` | the same at 2 x 1.875 km, run non-dimensional | `nondimensionalization = True` |
| `continental_rift_non_dimensional/model2`, `model3` | 3D in Mandyoc: the xz section | |
| `variable_velocity_field` | extension, convergence from 20 Myr and stop at 70 Myr (`scale_bcv.txt`) | accordion with `ti_convergence`, `ti_quiescence`, `quiescence_scale = 0` |
| `Crameri2012_case2` | plume rising below a lid, constant viscosity | isothermal (`temperature_difference = 0`) |
| `slab_detachment` | lithosphere with a vertical slab (Schmalholz, 2011) | `sticky_air = off`, isothermal |
| `vanKeken1997_case1a` | Rayleigh-Taylor instability (van Keken et al., 1997) | `sticky_air = off`, isothermal |
| `two_layers_with_variable_conductivity` | two layers with different thermal diffusivity | `sticky_air = off`, `preset_initial_temperature` (field) |
| `punch` | rigid punch (Prandtl), 1 x 0.5 m | `sticky_air = off`, `sampling = 0.001 m`, `custom_veloc_field = build_veloc.py` |

```bash
cd continental_rift
tapioca-svg continental_rift.svg -o inputs
```

Settings of the SVG scenarios used here (see `SVGScenario`):

- `ti_convergence`, `ti_quiescence` and `quiescence_scale` set the times of the accordion and the velocity of its
  quiescence (`scale_bcv.txt`). A `time_max` at the end of the quiescence has no second extension, and
  `quiescence_scale = 0` stops the boundaries.
- `sticky_air = off`: the top unit drawn is a solid. tapioca adds a sticky air of zero thickness at the top of the
  domain, as in the punch of Mandyoc, and writes `sticky_blanket_air = False`.
- `sampling`: the vertical sampling of the SVG that locates the interfaces (100 m by default), smaller than the
  thinnest unit and than the domain (1 mm in the punch).
- Any parameter of `param.txt` (keys of `tapioca.PARAMS_TYPES`) in the description of the setup shape replaces the
  value of tapioca. With `temperature_difference = 0` the scenario is isothermal: the initial temperature is 0 oC.
- `custom_veloc_field`: `punch/build_veloc.py` builds the velocity of the punch (`build_velocity`) and writes its
  `bcv_0.txt` (`write_files`, with `bcv_extern = True`).

What still differs from the notebooks:

- The initial temperature of the rifts is the one of tapioca (conduction with the LAB at 1350 oC). The notebooks
  use 1300 oC at the LAB, Tp = 1262 oC and 500 Myr of diffusion.
- `continental_rift_with_surface_processes`: the current Mandyoc has no `sp_mode = 5` (erosion with
  `precipitation.txt`), so it uses diffusion. The Mandyoc of ggciag reads `sediment_layers.txt` whenever
  `sp_surface_processes` is True, so tapioca adds a layer of sediments every 10 Myr.
- The sticky air of zero thickness added by tapioca (scenarios without sticky air) has the air material of tapioca.
  The second interface of the punch notebook has the material of the punch. Neither layer has particles.
- The figures of the punch show its lengths in km (0 km), since the domain is 1 m.
