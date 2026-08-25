import sys
import h5py
import numpy as np

if len(sys.argv) != 2:
    print("Usage: python makeAREPOsnap.py <3D_file.hdf5>")
    sys.exit(1)

InputFileName = sys.argv[1]

OutputFileName = InputFileName.replace(".hdf5", "_AREPO.hdf5")

gamma = 5.0 / 3.0
mu = 0.597
k_B = 1.380649e-16
m_p = 1.67262192369e-24

Unit_in_g = 1.989e43
Unit_in_cm = 3.085678e21
Unit_in_cm_per_s = 1.e5

Unit_in_erg = Unit_in_g * Unit_in_cm_per_s**2
Unit_in_erg_per_cm3 = Unit_in_erg / Unit_in_cm**3
Unit_in_G = np.sqrt(Unit_in_erg_per_cm3)
Unit_in_muG = Unit_in_G / 1e-6

Unit_density = Unit_in_g / Unit_in_cm**3

with h5py.File(InputFileName, "r") as f:
    Coordinates = f["Coordinates"][:]
    Density_cgs = f["Density"][:]
    Temperature = f["Temperature"][:]
    Velocity_cgs = f["Velocity"][:]
    Metallicity = f["Metallicity"][:]
    Volume_cgs = f["Volume"][:]
    MagneticField_cgs = f["MagneticField"][:] if "MagneticField" in f else None

    boxsize = f.attrs["boxsize"]/ Unit_in_cm
    redshift = f.attrs["redshift"]

NPart = len(Density_cgs)

Coordinates = Coordinates / Unit_in_cm + boxsize/2
Velocity = Velocity_cgs / Unit_in_cm_per_s
Density = Density_cgs / Unit_density
Volume = Volume_cgs / Unit_in_cm**3

mask = (
    (Coordinates[:, 0] >= 0) & (Coordinates[:, 0] <= 2*boxsize) &
    (Coordinates[:, 1] >= 0) & (Coordinates[:, 1] <= 2*boxsize) &
    (Coordinates[:, 2] >= 0) & (Coordinates[:, 2] <= 2*boxsize)
)

Coordinates = Coordinates[mask]
Velocity = Velocity[mask]
Density = Density[mask]
Volume = Volume[mask]

Mass_cgs = Density_cgs * Volume_cgs
Masses = Mass_cgs / (Unit_in_g)

## make PDFs
mass_bins = np.logspace(
    np.log10(Masses.min()),
    np.log10(Masses.max()),
    50
)

volume_bins = np.logspace(
    np.log10(Volume.min()),
    np.log10(Volume.max()),
    50
)

import matplotlib.pyplot as plt 

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Mass PDF
axes[0].hist(
    Masses,
    bins=mass_bins,
    density=True
)

axes[0].set_xscale("log")
axes[0].set_yscale("log")
axes[0].set_xlabel(r"Cell mass [$10^{10}~M_\odot$]")
axes[0].set_ylabel("PDF")
axes[0].set_title("Cell Mass PDF")

# Volume PDF
axes[1].hist(
    Volume,
    bins=volume_bins,
    density=True
)

axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].set_xlabel(r"Cell volume [$\mathrm{kpc}^3$]")
axes[1].set_ylabel("PDF")
axes[1].set_title("Cell Volume PDF")

plt.savefig(InputFileName.replace(".hdf5", "_MV_PDF.pdf"))


Utherm_cgs = k_B * Temperature / ((gamma - 1.0) * mu * m_p)
Utherm = Utherm_cgs / Unit_in_cm_per_s**2

if MagneticField_cgs is not None:
    MagneticField = MagneticField_cgs / Unit_in_G
else:
    MagneticField = None

with h5py.File(OutputFileName, "w") as IC:

    header = IC.create_group("Header")
    part0 = IC.create_group("PartType0")

    NumPart = np.array([NPart, 0, 0, 0, 0, 0], dtype=np.int32)

    header.attrs.create("NumPart_ThisFile", NumPart)
    header.attrs.create("NumPart_Total", NumPart)
    header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=np.int32))
    header.attrs.create("MassTable", np.zeros(6, dtype=np.float64))
    header.attrs.create("Time", 1.0 / (1.0 + redshift))
    header.attrs.create("Redshift", redshift)
    header.attrs.create("BoxSize", boxsize)
    header.attrs.create("NumFilesPerSnapshot", 1)
    header.attrs.create("Omega0", 0.0)
    header.attrs.create("OmegaB", 0.0)
    header.attrs.create("OmegaLambda", 0.0)
    header.attrs.create("HubbleParam", 0.0)
    header.attrs.create("Flag_Sfr", 0)
    header.attrs.create("Flag_Cooling", 0)
    header.attrs.create("Flag_StellarAge", 0)
    header.attrs.create("Flag_Metals", 1)
    header.attrs.create("Flag_Feedback", 0)
    header.attrs.create("Flag_DoublePrecision", 1)

    header.attrs.create(
        "center",
        np.array([boxsize/2, boxsize/2, boxsize/2])
    )

    part0.create_dataset(
        "ParticleIDs",
        data=np.arange(1, NPart + 1, dtype=np.int64)
    )

    part0.create_dataset(
        "Coordinates",
        data=Coordinates
    )

    part0.create_dataset(
        "Masses",
        data=Masses
    )

    part0.create_dataset(
        "Metallicity",
        data=Metallicity[mask]/0.02
    )

    part0.create_dataset(
        "Velocities",
        data=Velocity
    )

    part0.create_dataset(
        "InternalEnergy",
        data=Utherm
    )

    part0.create_dataset(
        "Density",
        data=Density
    )

    if MagneticField is not None:
        part0.create_dataset(
            "MagneticField",
            data=MagneticField[mask,:]
        )

print(f"Saved AREPO snapshot: {OutputFileName}")
print(f"Number of particles: {NPart}")
print(f"BoxSize: {boxsize} kpc")