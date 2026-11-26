"""Figures of the tapIOca handbook (docs/handbook), made from the examples of the repository.

    python docs/handbook/scripts/make_figures.py WORK_FOLDER [STATE_FOLDER]

WORK_FOLDER receives the scenarios built from examples/svg_scenarios (tapioca-svg, 5 km). STATE_FOLDER is an optional
finished run (with lithology.nc, temperature.nc and interfaces.txt) for the figures of tapioca-svg-from-state.
The SVG drawings are exported with Inkscape, when it is installed.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
WORK = Path(sys.argv[1]).resolve(); WORK.mkdir(parents=True, exist_ok=True)
STATE = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else None
EXAMPLES = ["rifting", "compression", "accordion", "quiescence", "rifting_sedimentation", "rifting_with_basin"]
TAPIOCA_SVG = str(Path(sys.executable).parent / "tapioca-svg")      # the tools of the same environment


def run(*args):
    print(" ".join(map(str, args)))
    subprocess.run(list(map(str, args)), check=True, stdout=subprocess.DEVNULL)


def export_svg(svg, png, dpi=60):
    inkscape = shutil.which("inkscape")
    if inkscape is None:
        print("inkscape not found:", png, "not made"); return
    run(inkscape, svg, "--export-type=png", f"--export-filename={png}", f"--export-dpi={dpi}",
        "--export-background=white", "--export-background-opacity=1")


# ---- scenarios of examples/svg_scenarios: drawings, numerical setup and initial temperature ----
for name in EXAMPLES:
    svg = ROOT / "examples/svg_scenarios" / name / f"{name}.svg"
    folder = WORK / name
    if not (folder / "numerical_setup.png").exists():
        run(TAPIOCA_SVG, svg, "-o", folder, "--dx", 5, "--dz", 5, "--no-scripts", "--no-zip")
    shutil.copy(folder / "numerical_setup.png", OUT / f"setup_{name}.png")
    shutil.copy(folder / "initial_temperature.png", OUT / f"temperature_{name}.png")
    if name in ("rifting", "compression"):
        export_svg(svg, OUT / f"svg_{name}.png")

# ---- post-processing: MandyocScen and the deviatoric stress of examples/continental_rift ----
from tapioca.scenClasses import MandyocScen
import tapioca.post_processing  # noqa: F401  (registers the postproc accessor)
from cmcrameri import cm

scen = MandyocScen(ROOT / "examples/continental_rift", variables=["density", "temperature", "viscosity", "velocity", "strain_rate"],
                   load_lithology=True)
scen.DTree.postproc.DeviatoricStressTensor()
mesh, lith = scen.DTree["/mesh/original"].ds, scen.DTree["/mesh/upscaled"].ds
times = lith.time.values[[0, len(lith.time) // 2, -1]]
fig, axs = plt.subplots(2, 2, figsize=(7.4, 4.0), constrained_layout=True)
km = lambda a: a / 1e3
for ax, t in zip(axs.flat[:3], times):
    l = lith.lithology.sel(time=t)
    ax.pcolormesh(km(l.x), km(l.z), l, cmap=cm.batlowS, shading="auto", rasterized=True)
    ax.set_title(f"lithology, {t:.1f} Myr", loc="left", fontsize=9); ax.set_ylim(-160, 0)
# tau_J2 = (tau_xx^2 + tau_zz^2)/2 + tau_xz^2 is in MPa^2 (its attribute says MPa): sqrt(J2) is a stress in MPa
j2 = np.sqrt(mesh.tau_J2.sel(time=times[-1]).where(mesh.density.sel(time=times[-1]) > 200))   # rock only
im = axs[1, 1].pcolormesh(km(j2.x), km(j2.z), j2, cmap=cm.lajolla, shading="auto", rasterized=True,
                          vmin=0, vmax=float(np.nanpercentile(j2, 99)))
axs[1, 1].set_title(f"sqrt(J2) of the deviatoric stress (MPa), {times[-1]:.1f} Myr", loc="left", fontsize=9); axs[1, 1].set_ylim(-160, 0)
plt.colorbar(im, ax=axs[1, 1], fraction=0.05)
for ax in axs.flat:
    ax.set_xlabel("x (km)", fontsize=8); ax.set_ylabel("z (km)", fontsize=8); ax.tick_params(labelsize=7)
fig.savefig(OUT / "postproc_continental_rift.png", dpi=200); plt.close(fig)
print(OUT / "postproc_continental_rift.png")

# ---- tapioca-svg-from-state: SVG of the final state of a run, and the new scenario ----
if STATE is not None:
    svg = WORK / "final_state.svg"
    run(Path(TAPIOCA_SVG).with_name("tapioca-svg-from-state"), STATE, "-o", svg)
    export_svg(svg, OUT / "svg_final_state.png", dpi=40)
    folder = WORK / "from_state"
    if not (folder / "numerical_setup.png").exists():
        run(TAPIOCA_SVG, svg, "-o", folder, "--no-scripts", "--no-zip")
    shutil.copy(folder / "numerical_setup.png", OUT / "setup_from_state.png")
    shutil.copy(folder / "initial_temperature.png", OUT / "temperature_from_state.png")
print("done")
