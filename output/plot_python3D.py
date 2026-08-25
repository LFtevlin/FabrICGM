import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy.units as un
import astropy.constants as const
from matplotlib.colors import LogNorm, SymLogNorm

if len(sys.argv) < 2:
    print("Usage: python plot_solution.py solution_file.hdf5")
    sys.exit()

filename = sys.argv[1]

with h5py.File(filename, "r") as f:
    coordinates = (f["Coordinates"][:] * un.cm).to(un.kpc)
    rho = f["Density"][:] * un.g / un.cm**3
    P = f["Pressure"][:] * un.erg / un.cm**3
    velocity = (f["Velocity"][:] * un.cm / un.s).to(un.km / un.s)
    B = f["MagneticField"][:] * un.G if "MagneticField" in f else None

x = coordinates[:, 0]
y = coordinates[:, 1]
z = coordinates[:, 2]

N = len(x)

print(f"Number of particles: {N}")
print("x:", x.min(), x.max())
print("y:", y.min(), y.max())
print("z:", z.min(), z.max())

vx = velocity[:, 0]
vy = velocity[:, 1]
vz = velocity[:, 2]

v = np.sqrt(vx**2 + vy**2 + vz**2)

r = np.sqrt(x**2 + y**2 + z**2)

vr = (x * vx + y * vy + z * vz) / np.maximum(r, 1e-30 * un.kpc)

if B is not None:
    Bx = B[:, 0]
    By = B[:, 1]
    Bz = B[:, 2]
    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
else:
    Bmag = None

gamma = 5.0 / 3.0
mu = 0.597

T = ((gamma - 1.0) * mu * const.m_p / const.k_B * P / rho).to(un.K)
T[rho == 0] = 0 * un.K

mask_z = np.abs(z.value) < 1
mask_y = np.abs(y.value) < 1

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
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
    "figure.dpi": 150,
})

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))

ax = axes[0]

ax.scatter(x[mask_z], y[mask_z], s=0.1, alpha=0.4, rasterized=True)

ax.set_xlabel(r"$x\ [{\rm kpc}]$")
ax.set_ylabel(r"$y\ [{\rm kpc}]$")
ax.set_aspect("equal", adjustable="box")

ax = axes[1]

ax.scatter(x[mask_y], z[mask_y], s=0.1, alpha=0.4, rasterized=True)

ax.set_xlabel(r"$x\ [{\rm kpc}]$")
ax.set_ylabel(r"$z\ [{\rm kpc}]$")
ax.set_aspect("equal", adjustable="box")

fig.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.97, wspace=0.25)

out = filename.replace(".hdf5", "_sampling.png")

plt.savefig(out, dpi=400, bbox_inches="tight")


N_bins = 50

r_bins = np.linspace(r.min(), r.max(), N_bins + 1)
r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

def radial_median(r, value, r_bins):
    result = np.zeros(len(r_bins) - 1) * value.unit

    for i in range(len(r_bins) - 1):
        mask = (r >= r_bins[i]) & (r < r_bins[i + 1])

        if np.any(mask):
            result[i] = np.median(value[mask])

    return result

rho_radial = radial_median(r, rho, r_bins)
T_radial = radial_median(r, T, r_bins)
v_radial = radial_median(r, v, r_bins)
vr_radial = radial_median(r, vr, r_bins)

if Bmag is not None:
    B_radial = radial_median(r, Bmag, r_bins)

Mdot = (4.0 * np.pi * r**2 * rho * vr).to(un.M_sun / un.yr)
Mdot_radial = radial_median(r, Mdot, r_bins)

fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.2))
axes = axes.flatten()

ax = axes[0]

ax.plot(r_centers, rho_radial)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\rho\ [{\rm g\,cm^{-3}}]$")

ax = axes[1]

ax.plot(r_centers, T_radial)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$T\ [{\rm K}]$")

ax = axes[2]

ax.plot(r_centers, v_radial)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$|v|\ [{\rm km\,s^{-1}}]$")

ax = axes[3]

ax.plot(r_centers, vr_radial)
ax.axhline(0, linestyle="--", linewidth=0.8)
ax.set_xscale("log")
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$v_r\ [{\rm km\,s^{-1}}]$")

ax = axes[4]

if Bmag is not None:
    ax.plot(r_centers, B_radial)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$|B|\ [{\rm G}]$")
else:
    ax.text(0.5, 0.5, "No magnetic field", ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel(r"$r$")

ax = axes[5]

ax.plot(r_centers, Mdot_radial)
ax.axhline(0, linestyle="--", linewidth=0.8)
ax.set_xscale("log")
ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$\dot{M}\ [{\rm M_\odot\,yr^{-1}}]$")

for ax in axes:
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.35)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.15)

fig.subplots_adjust(left=0.08, right=0.98, bottom=0.09, top=0.97, wspace=0.35, hspace=0.35)

out = filename.replace(".hdf5", "_profiles.png")

fig.savefig(out, dpi=400, bbox_inches="tight")

slice_thickness = 300 * un.kpc
npix = 100
vector_npix = 15
vector_slice = 5.0

cmaps = {
    "Density": "magma",
    "Temperature": "plasma",
    "Velocity": "RdBu_r",
    "Magnetic field": "viridis",
}

quantities = [
    ("Density", rho, r"$\rho_{\mathrm{turb}}\ [{\rm g\,cm^{-3}}]$"),
    ("Temperature", T, r"$T_{\mathrm{turb}}\ [{\rm K}]$"),
    ("Velocity", v, r"$|v_{\mathrm{turb}}|\ [{\rm km\,s^{-1}}]$"),
]

if Bmag is not None:
    quantities.append(("Magnetic field", Bmag, r"$|B_{\mathrm{turb}}|\ [{\rm G}]$"))

ncols = len(quantities)

def make_projection(xdata, ydata, los, value, weight, extent, los_halfwidth, npix=128):
    mask = np.isfinite(xdata) & np.isfinite(ydata) & np.isfinite(los) & np.isfinite(value) & np.isfinite(weight) & (weight > 0) & (value > 0) & (np.abs(los) < los_halfwidth)
    xdata = xdata[mask]
    ydata = ydata[mask]
    value = value[mask]
    weight = weight[mask]

    xmin, xmax, ymin, ymax = extent
    x_edges = np.linspace(xmin, xmax, npix + 1)
    y_edges = np.linspace(ymin, ymax, npix + 1)

    weighted_value, _, _ = np.histogram2d(xdata, ydata, bins=[x_edges, y_edges], weights=value * weight)
    total_weight, _, _ = np.histogram2d(xdata, ydata, bins=[x_edges, y_edges], weights=weight)
    projection = np.divide(weighted_value, total_weight, out=np.full_like(weighted_value, np.nan, dtype=float), where=total_weight > 0)

    return projection, x_edges, y_edges

xplot = x.to(un.kpc).value
yplot = y.to(un.kpc).value
zplot = z.to(un.kpc).value

rmax = max(np.nanmax(np.abs(xplot)), np.nanmax(np.abs(yplot)), np.nanmax(np.abs(zplot)))

extent = (-rmax / 2, rmax / 2, -rmax / 2, rmax / 2)

fig, axes = plt.subplots(2, ncols, figsize=(2.2 * ncols, 4.5), squeeze=False)

slice_kpc = slice_thickness.to(un.kpc).value

for col, (name, value, label) in enumerate(quantities):
    vv_all = np.asarray(value.value if hasattr(value, "value") else value, dtype=float)
    weight = np.ones_like(vv_all)

    proj_xy, x_edges, y_edges = make_projection(xplot, yplot, zplot, vv_all, weight, extent, slice_kpc, npix=npix)
    ax = axes[0, col]

    finite = np.isfinite(proj_xy) & (proj_xy > 0)

    if np.any(finite):
        vmin = np.nanpercentile(proj_xy[finite], 1)
        vmax = np.nanpercentile(proj_xy[finite], 99)
        im = ax.pcolormesh(x_edges, y_edges, proj_xy.T, shading="auto", cmap=cmaps[name], norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
        cbar = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.046)
        cbar.set_label(label)
        cbar.ax.tick_params(direction="in", labelsize=7, width=0.6, length=3)

    if col == 0:
        ax.set_ylabel(r"$y\ [{\rm kpc}]$")

    ax.set_aspect("equal", adjustable="box")

    proj_xz, x_edges, z_edges = make_projection(xplot, zplot, yplot, vv_all, weight, extent, slice_kpc, npix=npix)
    ax = axes[1, col]

    finite = np.isfinite(proj_xz) & (proj_xz > 0)

    if np.any(finite):
        vmin = np.nanpercentile(proj_xz[finite], 1)
        vmax = np.nanpercentile(proj_xz[finite], 99)
        im = ax.pcolormesh(x_edges, z_edges, proj_xz.T, shading="auto", cmap=cmaps[name], norm=LogNorm(vmin=vmin, vmax=vmax), rasterized=True)
        cbar = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.046)
        cbar.set_label(label)
        cbar.ax.tick_params(direction="in", labelsize=7, width=0.6, length=3)

    ax.set_xlabel(r"$x\ [{\rm kpc}]$")

    if col == 0:
        ax.set_ylabel(r"$z\ [{\rm kpc}]$")

    ax.set_aspect("equal", adjustable="box")

vxx = vx.to(un.km / un.s).value
vyy = vy.to(un.km / un.s).value
vzz = vz.to(un.km / un.s).value

if B is not None:
    Bxx = Bx.to(un.G).value
    Byy = By.to(un.G).value
    Bzz = Bz.to(un.G).value

def bin_vector_field(xdata, ydata, udata, vdata, extent, npix=18):
    xmin, xmax, ymin, ymax = extent
    x_edges = np.linspace(xmin, xmax, npix + 1)
    y_edges = np.linspace(ymin, ymax, npix + 1)

    mask = np.isfinite(xdata) & np.isfinite(ydata) & np.isfinite(udata) & np.isfinite(vdata)
    xdata = xdata[mask]
    ydata = ydata[mask]
    udata = udata[mask]
    vdata = vdata[mask]

    sum_u, _, _ = np.histogram2d(xdata, ydata, bins=[x_edges, y_edges], weights=udata)
    sum_v, _, _ = np.histogram2d(xdata, ydata, bins=[x_edges, y_edges], weights=vdata)
    counts, _, _ = np.histogram2d(xdata, ydata, bins=[x_edges, y_edges])

    mean_u = np.divide(sum_u, counts, out=np.full_like(sum_u, np.nan), where=counts > 0)
    mean_v = np.divide(sum_v, counts, out=np.full_like(sum_v, np.nan), where=counts > 0)

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    X, Y = np.meshgrid(x_centers, y_centers, indexing="ij")

    return X, Y, mean_u, mean_v, counts

def add_vectors(ax, xdata, ydata, udata, vdata, extent, npix=18, arrow_length=0.04):
    X, Y, U, V, counts = bin_vector_field(xdata, ydata, udata, vdata, extent, npix=npix)

    magnitude = np.sqrt(U**2 + V**2)

    valid = np.isfinite(U) & np.isfinite(V) & np.isfinite(magnitude) & (magnitude > 0) & (counts > 0)

    Xq = X[valid]
    Yq = Y[valid]
    Uq = U[valid]
    Vq = V[valid]
    magnitude = magnitude[valid]

    Udir = Uq / magnitude
    Vdir = Vq / magnitude

    width = extent[1] - extent[0]
    length = arrow_length * width

    Uplot = Udir * length
    Vplot = Vdir * length

    ax.quiver(Xq, Yq, Uplot, Vplot, color="black", angles="xy", scale_units="xy", scale=1, width=0.0035, headwidth=4.0, headlength=5.0, headaxislength=4.5, pivot="mid", rasterized=True)

mask_xy = np.abs(zplot) < vector_slice
mask_xz = np.abs(yplot) < vector_slice

velocity_col = 2

add_vectors(axes[0, velocity_col], xplot[mask_xy], yplot[mask_xy], vxx[mask_xy], vyy[mask_xy], extent, npix=vector_npix, arrow_length=0.055)
add_vectors(axes[1, velocity_col], xplot[mask_xz], zplot[mask_xz], vxx[mask_xz], vzz[mask_xz], extent, npix=vector_npix, arrow_length=0.055)

if B is not None:
    B_col = 3
    add_vectors(axes[0, B_col], xplot[mask_xy], yplot[mask_xy], Bxx[mask_xy], Byy[mask_xy], extent, npix=vector_npix, arrow_length=0.055)
    add_vectors(axes[1, B_col], xplot[mask_xz], zplot[mask_xz], Bxx[mask_xz], Bzz[mask_xz], extent, npix=vector_npix, arrow_length=0.055)

for ax in axes.flat:
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.25)

if B is None:
    axes[0, -1].axis("off")
    axes[1, -1].axis("off")

fig.subplots_adjust(left=0.07, right=0.98, bottom=0.08, top=0.94, wspace=0.65, hspace=0.15)

out = filename.replace(".hdf5", "_combined_projection.png")

fig.savefig(out, dpi=400, bbox_inches="tight")
plt.close(fig)


quantities = [
    ("vx", vx, r"$v_x\ [{\rm km/s}]$"),
    ("vy", vy, r"$v_y\ [{\rm km/s}]$"),
    ("vz", vz, r"$v_z\ [{\rm km/s}]$")
]

ncols = len(quantities)

fig, axes = plt.subplots(2, ncols, figsize=(2.5 * ncols, 4.8), squeeze=False)

slice_kpc = slice_thickness.to(un.kpc).value

for col, (name, value, label) in enumerate(quantities):
    if hasattr(value, "value"):
        vv_all = value.value
    else:
        vv_all = np.asarray(value)

    vv_all = np.asarray(vv_all, dtype=float)

    mask_xy = (np.abs(z) < slice_thickness)

    xx = np.asarray(x[mask_xy])
    yy = np.asarray(y[mask_xy])
    vv = vv_all[mask_xy]

    valid = np.isfinite(xx) & np.isfinite(yy) & np.isfinite(vv)

    xx = xx[valid]
    yy = yy[valid]
    vv = vv[valid]

    ax = axes[0, col]

    vmin = np.nanmin(vv)
    vmax = np.nanmax(vv)

    sc = ax.scatter(xx, yy, c=vv, s=0.5, alpha=0.7, cmap="RdBu_r", norm=SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax), rasterized=True)

    cbar = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.046)
    cbar.set_label(label)
    cbar.ax.tick_params(direction="in", labelsize=7, width=0.6, length=3)

    ax.set_xlabel(r"$x\ [{\rm kpc}]$")
    ax.set_ylabel(r"$y\ [{\rm kpc}]$")
    ax.set_title(rf"{name}: $x$-$y$" "\n" rf"$|z| < {slice_kpc:g}\,\mathrm{{kpc}}$", pad=4)
    ax.set_aspect("equal", adjustable="box")

    mask_xz = (np.abs(y) < slice_thickness)

    xx = np.asarray(x[mask_xz])
    zz = np.asarray(z[mask_xz])
    vv = vv_all[mask_xz]

    valid = np.isfinite(xx) & np.isfinite(zz) & np.isfinite(vv)

    xx = xx[valid]
    zz = zz[valid]
    vv = vv[valid]

    ax = axes[1, col]

    vmin = np.nanmin(vv)
    vmax = np.nanmax(vv)

    sc = ax.scatter(xx, zz, c=vv, s=0.5, alpha=0.7, cmap="RdBu_r", norm=SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax), rasterized=True)

    cbar = fig.colorbar(sc, ax=ax, pad=0.015, fraction=0.046)
    cbar.set_label(label)
    cbar.ax.tick_params(direction="in", labelsize=7, width=0.6, length=3)

    ax.set_xlabel(r"$x\ [{\rm kpc}]$")
    ax.set_ylabel(r"$z\ [{\rm kpc}]$")
    ax.set_title(rf"{name}: $x$-$z$" "\n" rf"$|y| < {slice_kpc:g}\,\mathrm{{kpc}}$", pad=4)
    ax.set_aspect("equal", adjustable="box")

fig.subplots_adjust(left=0.07, right=0.98, bottom=0.08, top=0.94, wspace=0.55, hspace=0.35)

out = filename.replace(".hdf5", "_scatter_projection_velocities.png")

fig.savefig(out, dpi=400, bbox_inches="tight")

plt.close(fig)


R = np.sqrt(x**2 + y**2)

vtheta = ((x * vx + y * vy) * z - (x**2 + y**2) * vz) / (r * R)

vphi = (-y * vx + x * vy) / R

components = [
    (vr, r"$v_r$", "tab:blue"),
    (vphi, r"$v_\phi$", "tab:orange"),
    (vtheta, r"$v_\theta$", "tab:green")
]

fig, ax = plt.subplots(figsize=(3.6, 3.2))

for data, label, color in components:
    if hasattr(data, "value"):
        data = data.value

    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]

    pdf, bins = np.histogram(data, bins=50, density=True)
    centers = 0.5 * (bins[:-1] + bins[1:])

    ax.plot(centers, pdf, lw=1.3, color=color, label=label)

ax.set_xlabel(r"Velocity $[{\rm km\,s^{-1}}]$")
ax.set_ylabel(r"PDF")

ax.set_yscale("log")

ax.legend(frameon=False, loc="best", handlelength=2.0)

ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.35)
ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.15)

fig.subplots_adjust(left=0.16, right=0.97, bottom=0.16, top=0.96)

out = filename.replace(".hdf5", "_velocity_PDFs.png")

fig.savefig(out, dpi=400, bbox_inches="tight")

plt.close(fig)

print("Saved 3D pictures.")
