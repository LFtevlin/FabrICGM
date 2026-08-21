import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print("Usage: python plot_solution.py solution_file.hdf5")
    sys.exit(1)

filename = sys.argv[1]

kpc = 3.085677581e21
Msun = 1.98847e33

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "lines.linewidth": 1.2,
    "figure.dpi": 150,
})

with h5py.File(filename, "r") as f:
    radius = np.asarray(f["Radius"][:], dtype=float)
    density = np.asarray(f["Density"][:], dtype=float)
    temperature = np.asarray(f["Temperature"][:], dtype=float)
    velocity = np.asarray(f["Velocity"][:], dtype=float)
    metallicity = np.asarray(f["Metallicity"][:], dtype=float)
    Lambda = np.asarray(f["Lambda"][:], dtype=float)
    mass = np.asarray(f["Mass"][:], dtype=float)
    potential = np.asarray(f["Potential"][:], dtype=float)
    attrs = dict(f.attrs)

r_kpc = radius / kpc
mass_msun = mass / Msun
M200_msun = attrs["M200"] / Msun
Mgas_msun = attrs["Mgas"] / Msun
Mstar_msun = attrs["Mstellar"] / Msun
potential_km2s2 = potential / 1e10
Rsonic = attrs["Rsonic"] / kpc

if "R200" in attrs:
    R200 = attrs["R200"] / kpc
else:
    print("Warning: R200 not found in attributes")
    R200 = None

if R200 is not None:
    r_normalized = r_kpc / R200
else:
    r_normalized = r_kpc

if "Mcgm_R200" in attrs:
    Mcgm = attrs["Mcgm_R200"] / Msun
else:
    Mcgm = np.nan

if "Mdot" in attrs:
    Mdot = attrs["Mdot"] / Msun * 365.25 * 24.0 * 3600.0
else:
    Mdot = np.nan

fig, axes = plt.subplots(4, 2, figsize=(7.2, 8.0), sharex="col")
axes = axes.flatten()

def mark_radii(ax):
    if R200 is not None:
        Rsonic_normalized = Rsonic / R200
        ax.axvline(Rsonic_normalized, color="tab:red", linestyle="--", linewidth=0.9, label=r"$R_{\rm sonic}$", zorder=1)
        ax.axvline(1.0, color="black", linestyle=":", linewidth=0.9, label=r"$R_{200}$", zorder=1)
    else:
        ax.axvline(Rsonic, color="tab:red", linestyle="--", linewidth=0.9, label=r"$R_{\rm sonic}$", zorder=1)

ax = axes[0]
ax.plot(r_normalized, temperature, linewidth=1.3)
ax.set_ylabel(r"$T\ [{\rm K}]$")
ax.set_yscale("log")
mark_radii(ax)

ax = axes[1]
ax.plot(r_normalized, density, linewidth=1.3)
ax.set_ylabel(r"$\rho\ [{\rm g\,cm^{-3}}]$")
ax.set_yscale("log")
mark_radii(ax)

ax = axes[2]
ax.plot(r_normalized, velocity / 1e5, linewidth=1.3)
ax.set_ylabel(r"$v\ [{\rm km\,s^{-1}}]$")
mark_radii(ax)

ax = axes[3]
ax.plot(r_normalized, metallicity, linewidth=1.3)
ax.set_ylabel(r"$Z/Z_\odot$")
mark_radii(ax)

ax = axes[4]
ax.plot(r_normalized, Lambda, linewidth=1.3)
ax.set_ylabel(r"$\Lambda\ [{\rm erg\,cm^3\,s^{-1}}]$")
ax.set_yscale("log")
mark_radii(ax)

ax = axes[5]
ax.plot(r_normalized, mass_msun, linewidth=1.3)
ax.set_ylabel(r"$M(<r)\ [M_\odot]$")
ax.set_yscale("log")
mark_radii(ax)

ax = axes[6]
ax.plot(r_normalized, potential_km2s2, linewidth=1.3)
ax.set_ylabel(r"$\Phi\ [{\rm km^2\,s^{-2}}]$")
mark_radii(ax)

ax = axes[7]
ax.axis("off")

text_lines = [
    rf"$M_{{200}} = {M200_msun:.3e}\ M_\odot$",
    rf"$M_{{\rm gas}} = {Mgas_msun:.3e}\ M_\odot$",
    rf"$M_\star = {Mstar_msun:.3e}\ M_\odot$",
    rf"$M_{{\rm CGM}}(<R_{{200}}) = {Mcgm:.3e}\ M_\odot$",
]

if "fbaryon_R200" in attrs:
    text_lines.append(rf"$f_{{\rm B}}(<R_{{200}}) = {attrs['fbaryon_R200']:.3f}$")

text_lines.extend([
    rf"$\dot{{M}} = {Mdot:.3e}\ M_\odot\,\mathrm{{yr}}^{{-1}}$",
    rf"$R_{{\rm sonic}} = {Rsonic:.2f}\ \mathrm{{kpc}}$",
])

if R200 is not None:
    text_lines.append(rf"$R_{{200}} = {R200:.2f}\ \mathrm{{kpc}}$")

text = "\n\n".join(text_lines)

ax.text(0.05, 0.95, text, transform=ax.transAxes, va="top", ha="left", fontsize=8, linespacing=1.3)

axes[6].set_xlabel(r"$r/R_{200}$" if R200 is not None else r"$r\ [{\rm kpc}]$")

for ax in axes[:7]:
    ax.set_xscale("log")

    if R200 is not None:
        ax.set_xlim(5e-2, 3.7)

    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.35)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.15)

    handles, labels = ax.get_legend_handles_labels()

    if handles:
        ax.legend(frameon=False, loc="best", handlelength=2.0)

fig.subplots_adjust(left=0.10, right=0.98, bottom=0.07, top=0.98, wspace=0.28, hspace=0.22)

out = filename.replace(".hdf5", "_profiles.png")

fig.savefig(out, dpi=400, bbox_inches="tight")

print(f"Saved {out}")