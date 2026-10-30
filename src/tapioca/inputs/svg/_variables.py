'''
Variables of the SVG scenarios (tapioca.inputs.svg.SVGScenario): regimes, parameters of param.txt and the
scripts written with the input files.
'''

import os

__all__ = ["REGIMES", "MANTLE_RHEOLOGIES", "SYSTEMS", "DEFAULT_PARAMS",
           "ORGANIZE_SCRIPT", "RUN_SCRIPT", "POST_PROCESSING_SCRIPT", "NO_POST_PROCESSING_SCRIPT",
           "RUN_LAUNCHER", "ZIP_SCRIPT"]

REGIMES:list = ['rifting', 'compression', 'accordion', 'quiescence']
"""list: Kinematic regimes of the SVG scenarios."""
MANTLE_RHEOLOGIES:list = ['wet_olivine', 'dry_olivine']
"""list: Rheologies of the lithospheric mantle (its thickness is the reference of the velocity with h_v_reference = 'lithospheric_mantle')."""
SYSTEMS:list = ['linux', 'mac']
"""list: Systems of the run scripts."""

#Script of tapioca that moves the outputs of Mandyoc to one folder per variable (copied to the scenario)
ORGANIZE_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'postrunning', '0_organize_outputs.sh')

#Run script of Mandyoc (@KEY@ are replaced in write_scripts). PETSC_DIR, PETSC_ARCH and MANDYOC
#set in the environment override the values of the script.
RUN_SCRIPT = '''#!/usr/bin/env bash
# Runs Mandyoc and moves the outputs to one folder per variable (0_organize_outputs.sh). Created by tapioca.
# PETSC_DIR, PETSC_ARCH and MANDYOC set in the environment override the values below.

PETSC_DIR=${PETSC_DIR:-@PETSC_DIR@}
PETSC_ARCH=${PETSC_ARCH-@PETSC_ARCH@}
MANDYOC=${MANDYOC:-@MANDYOC@}
NUMBER_OF_CORES=@NCORES@
MANDYOC_OPTIONS="@MANDYOC_OPTIONS@"

touch FD.out
@LAUNCHER@$PETSC_DIR/$PETSC_ARCH/bin/mpirun -n $NUMBER_OF_CORES $MANDYOC $MANDYOC_OPTIONS | tee FD.out

# Moves the outputs to one folder per variable
bash 0_organize_outputs.sh
@POST_PROCESSING@'''

#End of the run script without post-processing (post_processing=False)
NO_POST_PROCESSING_SCRIPT = '''
# Zips the inputs and the outputs (skips the files already in the zip)
# bash zipper.sh
'''

#Post-processing of the run script (@POST_PROCESSING@ of RUN_SCRIPT, NO_POST_PROCESSING_SCRIPT with
#post_processing=False): the .nc files with the Julia scripts of tapioca, the snapshots and the video of the
#scenario with plot_snapshots.py and frames_generator.py of post_processing (https://github.com/jopmas/post_processing),
#then the zip of the scenario (zipper.sh) and the removal of the output files that are in the zip.
POST_PROCESSING_SCRIPT = '''
# Post-processing: .nc files (Julia scripts of tapioca), snapshots and video of the scenario (plot_snapshots.py
# and frames_generator.py of post_processing), zip of the scenario and removal of the output files.
# TAPIOCA_DIR, POST_PROCESSING_DIR, JULIA and PYTHON set in the environment override the values below.
# A step is skipped, with a message, if what it needs is not found.
TAPIOCA_DIR=${TAPIOCA_DIR:-@TAPIOCA_DIR@}
POST_PROCESSING_DIR=${POST_PROCESSING_DIR:-@POST_PROCESSING_DIR@}
JULIA=${JULIA:-julia}
PYTHON=${PYTHON:-python}

if command -v "$JULIA" > /dev/null && [ -f "$TAPIOCA_DIR/src/postrunning/meshes2netcdf.jl" ]; then
    "$JULIA" -t $NUMBER_OF_CORES "$TAPIOCA_DIR/src/postrunning/meshes2netcdf.jl" "$PWD"
    if [ -d steps ]; then
        # last step as reference step: the tracking scripts of post_processing select the particles there
        LAST_STEP=$(find time -maxdepth 1 -name "time_*.txt" | sed -E 's|.*time_([0-9]+)\\.txt|\\1|' | grep -E '^[0-9]+$' | sort -n | tail -1)
        "$JULIA" -t $NUMBER_OF_CORES "$TAPIOCA_DIR/src/postrunning/particles2netcdf.jl" "$PWD" "$LAST_STEP"
    fi
else
    echo "Julia ($JULIA) or $TAPIOCA_DIR/src/postrunning not found: the .nc files were not created."
fi

if [ -f lithology.nc ] && [ -f "$POST_PROCESSING_DIR/plot_snapshots.py" ]; then
    "$PYTHON" "$POST_PROCESSING_DIR/plot_snapshots.py"
    "$PYTHON" "$POST_PROCESSING_DIR/frames_generator.py"
else
    echo "lithology.nc or $POST_PROCESSING_DIR not found: the snapshots and the video were not made."
fi

# Zips the inputs and the outputs (skips the files already in the zip)
bash zipper.sh

# Removes the output files (folders of 0_organize_outputs.sh and bc_velocity_*.txt) only if the zip is fine and
# has all of them. The inputs, the .nc files and the figures are kept. The lists are made with find and compared
# in files, never expanded in the command line, to avoid the "argument list too long" error.
ZIP="$(basename "$PWD").zip"
OUTPUT_FOLDERS="time density heat pressure surface strain_rate strain temperature viscosity velocity thermal_diffusivity X_depletion Phi dPhi lithos steps"
LISTS=$(mktemp -d)
for folder in $OUTPUT_FOLDERS; do
    if [ -d "$folder" ]; then find "$folder" -type f; fi
done > "$LISTS/outputs"
find . -maxdepth 1 -type f -name "bc_velocity_*.txt" | sed 's|^\\./||' >> "$LISTS/outputs"
LC_ALL=C sort -u -o "$LISTS/outputs" "$LISTS/outputs"

if [ -f "$ZIP" ] && unzip -tq "$ZIP" > /dev/null; then
    zipinfo -1 "$ZIP" | LC_ALL=C sort -u > "$LISTS/zipped"
    MISSING=$(LC_ALL=C comm -23 "$LISTS/outputs" "$LISTS/zipped" | wc -l | tr -d ' ')
    if [ "$MISSING" -eq 0 ]; then
        for folder in $OUTPUT_FOLDERS; do
            if [ -d "$folder" ]; then rm -rf "./$folder"; fi
        done
        find . -maxdepth 1 -type f -name "bc_velocity_*.txt" -delete
        echo "Output files removed: they are in $ZIP"
    else
        echo "$MISSING output file(s) are not in $ZIP: the output files were kept."
    fi
else
    echo "$ZIP was not created or is damaged: the output files were kept."
fi
rm -rf "$LISTS"
'''

#caffeinate keeps the Mac awake while Mandyoc runs
RUN_LAUNCHER = {'linux': '', 'mac': 'caffeinate -i '}

#Zips the scenario in <folder>.zip, skipping the files already in the zip. The outputs are in the
#folders created by 0_organize_outputs.sh.
ZIP_SCRIPT = '''#!/usr/bin/env bash
# Zips the files of Mandyoc in $DIRNAME.zip, skipping the files already in the zip. Created by tapioca.
set -eo pipefail

DIRNAME=${DIRNAME:-$(basename "$PWD")}
ZIP="${DIRNAME:?DIRNAME is not set}.zip"

# Patterns: files or folders (the folders are zipped with all their files)
patterns=(
"interfaces.txt"
"param.txt"
"input*_0.txt"
"vel_bc.txt"
"velz_bc.txt"
"scale_bcv.txt"
"infos_*.txt"
"*.svg"
"run*.sh"
"bc_velocity_*.txt"
"sp_surface_global_*.txt"
"FD.out"
"density"
"heat"
"pressure"
"surface"
"lithos"
"strain"
"strain_rate"
"temperature"
"time"
"velocity"
"viscosity"
"thermal_diffusivity"
"steps"
"Phi"
"dPhi"
"X_depletion"
"*.bin*.txt"
"bc*-1.txt"
"*.log"
# "_*.nc"
)

export LC_ALL=C   # same order for comm
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Files that match the patterns; in the folders, all their files
# (find avoids the "argument list too long" error)
for pat in "${patterns[@]}"; do
    find . -maxdepth 1 -name "$pat" ! -name "$ZIP" -exec find {} -type f \\;
done | sed 's|^\\./||' | sort -u > "$tmp/all"

# Files already in the zip
if [[ -f "$ZIP" ]]; then
    zipinfo -1 "$ZIP" | sort -u > "$tmp/zipped"
    grow=-g   # appends to the zip without rewriting it
else
    : > "$tmp/zipped"
    grow=
fi

# Only the files that are not in the zip yet
comm -23 "$tmp/all" "$tmp/zipped" > "$tmp/new"
n=$(wc -l < "$tmp/new" | tr -d ' ')

if [[ "$n" -eq 0 ]]; then
    echo "No new files for $ZIP"
else
    echo "Adding $n file(s) to $ZIP"
    zip -q $grow "$ZIP" -@ < "$tmp/new"
fi
'''

DEFAULT_PARAMS:dict = {
    # Simulation options
    "multigrid": 1,
    "solver": "direct",
    "denok": 1e-14,
    "particles_per_element": 40,
    "surface_particles_per_element": 40,
    "particles_perturb_factor": 0.7,
    "rtol": 1.0e-7,
    "RK4": "Euler",
    "Xi_min": 1.0e-5,
    "random_initial_strain": 0.3,
    "pressure_const": -1.0,
    "initial_dynamic_range": True,
    "periodic_boundary": False,
    "high_kappa_in_asthenosphere": False,
    "K_fluvial": 2.0e-7,
    "m_fluvial": 1.0,
    "sea_level": 0.0,
    "basal_heat": 0.0,
    # Surface processes
    "sp_surface_tracking": True,
    "sp_surface_processes": False,
    "plot_sediment": False,
    "a2l": True,
    "free_surface_stab": True,
    "theta_FSSA": 0.5,
    # Time constraints
    "step_max": 800000,
    "dt_max": 10000.0,
    "step_print": 50,
    "sub_division_time_step": 0.5,
    "initial_print_step": 0,
    "initial_print_max_time": 1.0e6,
    # Viscosity
    "viscosity_reference": 1.0e26,
    "viscosity_max": 1.0e25,
    "viscosity_min": 1.0e18,
    "viscosity_per_element": "constant",
    "viscosity_mean_method": "arithmetic",
    "viscosity_dependence": "pressure",
    # External ASCII inputs/outputs
    "interfaces_from_ascii": True,
    "temperature_from_ascii": True,
    "binary_output": False,
    "sticky_blanket_air": True,
    "print_step_files": True,
    "checkered": False,
    "geoq": "on",
    "geoq_fac": 100.0,
    # Physical parameters
    "temperature_difference": 1500.0,
    "thermal_expansion_coefficient": 3.28e-5,
    "thermal_diffusivity_coefficient": 1.0e-6,
    "gravity_acceleration": 10.0,
    "density_mantle": 3300.0,
    "external_heat": 0.0e-12,
    "heat_capacity": 1250.0,
    "non_linear_method": "on",
    "adiabatic_component": "on",
    "radiogenic_component": "on",
    "export_lithology": True,
    # Velocity boundary conditions
    "top_normal_velocity": "fixed",
    "top_tangential_velocity": "free",
    "bot_normal_velocity": "fixed",
    "bot_tangential_velocity": "free",
    "left_normal_velocity": "fixed",
    "left_tangential_velocity": "fixed",
    "right_normal_velocity": "fixed",
    "right_tangential_velocity": "fixed",
    "surface_velocity": 0.0e-2,
    "multi_velocity": False,
    # Temperature boundary conditions
    "top_temperature": "fixed",
    "bot_temperature": "fixed",
    "left_temperature": "fixed",
    "right_temperature": "fixed",
    "rheology_model": 19,
    "T_initial": 3,
}
"""dict: Parameters of param.txt of the SVG scenarios that do not depend on the scenario (any of them can be
changed with keyword arguments of SVGScenario.write_all or SVGScenario.parameters)."""
