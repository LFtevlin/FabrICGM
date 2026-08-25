import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm
from scipy.interpolate import RegularGridInterpolator
import astropy.units as un
import sys

if len(sys.argv) != 3:
    print("Usage: python plot.py NAME z")
    sys.exit(1)

NAME = sys.argv[1]
z = float(sys.argv[2])

print(f"NAME = {NAME}")
print(f"Redshift = {z}")

nbins = 200
xmin, xmax = -200, 200
ymin, ymax = -200, 200
thickness = 50.0
zmax = thickness / 2
xedges = np.linspace(xmin, xmax, nbins + 1)
yedges = np.linspace(ymin, ymax, nbins + 1)
extent = [xmin, xmax, ymin, ymax]

with h5py.File(f"./output/turbulence_{NAME}.hdf5", "r") as f:
    c1 = (f["Coordinates"][()] * un.cm).to(un.kpc).value
    P1 = f["Pressure"][()]

gamma = 5 / 3
u_th1 = P1 / (gamma - 1)

mask1 = ((c1[:, 0] >= xmin) & (c1[:, 0] <= xmax) &
         (c1[:, 1] >= ymin) & (c1[:, 1] <= ymax))

projection1, _, _ = np.histogram2d(c1[mask1, 0], c1[mask1, 1], bins=[xedges, yedges], weights=u_th1[mask1])

with h5py.File(f"./output/solution3D_{NAME}.hdf5", "r") as f:
    c2 = (f["Coordinates"][()] * un.cm).to(un.kpc).value
    P2 = f["Pressure"][()]

u_th2 = P2 / (gamma - 1)

mask2 = ((c2[:, 0] >= xmin) & (c2[:, 0] <= xmax) &
         (c2[:, 1] >= ymin) & (c2[:, 1] <= ymax))

projection2, _, _ = np.histogram2d(c2[mask2, 0], c2[mask2, 1], bins=[xedges, yedges], weights=u_th2[mask2])

with h5py.File(f"output/solution3D_{NAME}.hdf5", "r") as f:
    c = (f["Coordinates"][()] * un.cm).to(un.kpc).value
    T = f["Temperature"][()]
    rho = f["Density"][()]
    Z = f["Metallicity"][()]

with h5py.File("../cooling/UVB_dust1_CR1_G1_shield0.hdf5", "r") as f:
    cool = f["Tdep"]["Cooling"][()]
    heat = f["Tdep"]["Heating"][()]
    cool_ids = [x.decode().strip() for x in f["IdentifierCooling"][()]]
    heat_ids = [x.decode().strip() for x in f["IdentifierHeating"][()]]
    nHBins = 10**f["TableBins"]["DensityBins"][()]
    ZBins = 10**f["TableBins"]["MetallicityBins"][()]
    zBins = f["TableBins"]["RedshiftBins"][()]
    TBins = 10**f["TableBins"]["TemperatureBins"][()]

cool_prim = cool_ids.index("TotalPrim")
cool_met = cool_ids.index("TotalMetal")
heat_prim = heat_ids.index("TotalPrim")
heat_met = heat_ids.index("TotalMetal")

m_p = 1.67262192369e-24
X_H = 0.76
nH = X_H * rho / m_p

ind_z = np.argmin(np.abs(zBins - z))

cool_table = 10**cool[ind_z, :, :, :, cool_prim] + 10**cool[ind_z, :, :, :, cool_met]
heat_table = 10**heat[ind_z, :, :, :, heat_prim] + 10**heat[ind_z, :, :, :, heat_met]

cool_interp = RegularGridInterpolator((TBins, ZBins, nHBins), np.log10(cool_table), bounds_error=False, fill_value=np.nan)
heat_interp = RegularGridInterpolator((TBins, ZBins, nHBins), np.log10(heat_table), bounds_error=False, fill_value=np.nan)

points = np.column_stack([T, Z, nH])

cooling = 10**cool_interp(points)
heating = 10**heat_interp(points)
net = cooling - heating
cooling_rate = net * nH**2

mask3 = ((np.abs(c[:, 2]) < zmax) &
         (np.abs(c[:, 0]) < xmax) &
         (np.abs(c[:, 1]) < ymax) &
         np.isfinite(cooling_rate))

x_proj = c[mask3, 0]
y_proj = c[mask3, 1]
cool_proj = cooling_rate[mask3]

sum_cooling, _, _ = np.histogram2d(x_proj, y_proj, bins=[xedges, yedges], weights=cool_proj)
count, _, _ = np.histogram2d(x_proj, y_proj, bins=[xedges, yedges])

mean_cooling = np.divide(sum_cooling, count, out=np.full_like(sum_cooling, np.nan), where=count > 0)

fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)

finite1 = np.isfinite(projection1)
vmax1 = np.nanmax(np.abs(projection1[finite1]))

im1 = axes[0].imshow(projection1.T, origin="lower", extent=extent, aspect="equal", cmap="seismic",
                     norm=SymLogNorm(linthresh=1e-10, vmin=-vmax1, vmax=vmax1))
axes[0].set_xlabel("x [kpc]")
axes[0].set_ylabel("y [kpc]")
fig.colorbar(im1, orientation="horizontal", location="top", pad=0.08, ax=axes[0],
             label=r"Projected turbulent thermal energy density [erg/cm$^3$]")

finite2 = np.isfinite(projection2)
vmin2 = np.nanmin(projection2[finite2])
vmax2 = np.nanmax(projection2[finite2])

im2 = axes[1].imshow(projection2.T, origin="lower", extent=extent, aspect="equal", cmap="viridis",
                     norm=LogNorm(vmin=vmin2, vmax=vmax2))
axes[1].set_xlabel("x [kpc]")
axes[1].set_ylabel("y [kpc]")
fig.colorbar(im2, orientation="horizontal", location="top", pad=0.08, ax=axes[1],
             label=r"Projected turbulent energy density [erg/cm$^3$]")

finite3 = np.isfinite(mean_cooling)
vmax3 = np.nanmax(np.abs(mean_cooling[finite3]))

im3 = axes[2].imshow(mean_cooling.T, origin="lower", extent=extent, aspect="equal", cmap="seismic",
                     norm=SymLogNorm(linthresh=1e-32, vmin=-vmax3, vmax=vmax3))
axes[2].set_xlabel("x [kpc]")
axes[2].set_ylabel("y [kpc]")
fig.colorbar(im3, ax=axes[2], orientation="horizontal", location="top", pad=0.08,
             label=r"Mean volumetric cooling rate [erg/cm$^3$/s]")

plt.savefig(f"./output/turb_e_proj_{NAME}.pdf", bbox_inches="tight")
plt.close()