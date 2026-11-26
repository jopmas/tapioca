# Post running scripts

These scripts contain code to facilitate the data formatting that TapIOca expects. Most of these scripts were made to handle with .txt files. Note that .txt files will be become deprecated in future versions of Mandyoc.

The common routine is (1) to organise text file outputs and (2) to convert them into netcdf files.

An example:

```bash
bash $PATH_TO_TAPIOCA/tapioca/src/postrunning/0_organize_outputs.sh
julia -t THREADS $PATH_TO_TAPIOCA/tapioca/src/postrunning/meshes2netcdf.jl $SCENARIO_PATH
```

or using a slurm system:
```bash
bash $PATH_TO_TAPIOCA/tapioca/src/postrunning/0_organize_outputs.sh
srun -N 1 -c THREADS --mem=RAM julia -t THREADS $PATH_TO_TAPIOCA/tapioca/src/postrunning/meshes2netcdf.jl $SCENARIO_PATH
```

in which `THREADS` is the number of desired threads, `RAM` is required memory, and `SCENARIO_PATH` is the path to the modelled scenario folder.
To create the dataset files, you can filter variables in the sticky air layer by setting the variable `AIR_DENSITY_THRESHOLD` with the air density used (usually 1-10 kg/m³). This will set all densities lower than this threshold to zero.
The user can specify a number of chuncks to process the outputs using `CHUNKS`. It is useful for fine grids that can consume dozens of RAM memory. Additionaly, the user is free to set the data types for each variable using the dict `DTYPES`. 

The amount of necessary RAM can be estimated using the velocity and lithology datasets. Velocity is the most memory-intensive standard variable because it requires two 8 bytes arrays (vx and vz), and reading it also requires loading the density array (4 bytes) into memory simultaneously. 

- Buffer Memory: $\text{Nx} \times \text{Nz} \times \text{Steps per Chunk} \times 16 \text{ bytes}$
- Thread Overhead: $\text{Active Threads} \times \text{Nx} \times \text{Nz} \times 20 \text{ bytes}$ (Loading vx, vz, and density into the thread)
- Total: $\text{Max RAM (Velocity)} = \frac{\text{Nx} \times \text{Nz}}{1024^2} \times (16 \times \text{Steps per Chunk} + 20 \times \text{Active Threads})$ MB.

Lithology uses 1 byte by default (Int8) which is very light, but the output is upscaled by a factor of 5 in both axis. This means the lithology grid has roughly 25 times more elements than your original mesh. 
- Upscaled Elements: $\text{Nxl} \times \text{Nzl} \approx 25 \times (\text{Nx} \times \text{Nz})$
- Buffer Memory: $\text{Nxl} \times \text{Nzl} \times \text{Steps per Chunk} \times 1 \text{ byte}$
- Thread Overhead: $\text{Active Threads} \times \text{Nxl} \times \text{Nzl} \times 1 \text{ byte}$ (Plus a small CSV parsing overhead per thread)
- Total: $\text{Max RAM (Lithology)} = \frac{\text{Nxl} \times \text{Nzl}}{1024^2} \times (\text{Steps per Chunk} + \text{Active Threads})$ MB.

In which $\text{Steps per Chunk} = \text{Total Steps} / \text{CHUNKS}$. 

A practical amount of necessary RAM for each chunk is: $1.25 \times \max{R_{\text{litho}},R_{\text{veloc}}}$. Note that the increasing of the chunks also increases the processing time due to the writing data. For most cases (1201 x 301 elements with ~320 time steps), 2 chunks brings an optimized memory usage and a good time performance. 

## Lagrangian Particles (Markers)

The particles2netcdf.jl script converts Lagrangian markers into a NetCDF file. The output is structured with `(id, time)` dimensions to store `x`, `z`, and `lithology` variables, following the CF-Conventions for trajectory data (where each particle maintains a unique ID over time). The basic usage is:

```bash
srun -N 1 -c THREADS --mem=RAM julia -t THREADS $PATH_TO_TAPIOCA/tapioca/src/postrunning/particles2netcdf.jl $SCENARIO_PATH REFERENCE_TIME_STEP
```

The `REFERENCE_TIME_STEP` is the time step that will be used to get the required IDs. If this is the final step of the simulation, only the particles within the domain at the end will be tracked. On the other hand, if it is 0, all the initial particles will be tracked.

## Zarr outputs of Mandyoc

With `zarr_output = True` in `param.txt` (Mandyoc of the branch `zarr_outputs` of
[jopmas/mandyoc_jpm](https://github.com/jopmas/mandyoc_jpm)), Mandyoc writes one Zarr v2 store per variable, appended
at each printed step, in the layout of mandyocjll: `density.zarr`, ..., `velocity.zarr` (`vx`, `vz`),
`lithology.zarr`, `particles_trajectories.zarr` and `surface.zarr`, with the names, dimensions, units and types of
the `.nc` files, so they can be read directly (`xr.open_zarr("temperature.zarr")`), also while the run goes on.
There are no text outputs to organize. `zarr2netcdf.py` writes `<variable>.nc` from each store, with the compression
of the Julia scripts, and `times.csv`:

```bash
python $PATH_TO_TAPIOCA/tapioca/src/postrunning/zarr2netcdf.py $SCENARIO_PATH [variable ...]
```

It needs `zarr` and `h5netcdf` (environment `tapioca`); the times are read in parallel threads (`ZARR2NC_THREADS`,
default: the number of cores). The trajectories follow the particles printed in step 0 (`reference_timestep = 0`, as
in mandyocjll) and also have their `strain`. The Julia scripts are unchanged and keep converting the text outputs.
