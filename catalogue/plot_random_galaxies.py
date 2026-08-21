import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

G = 4.30091e-6
H0 = 70.0

def R_200(M200, z):
    Hz = H0 * np.sqrt(0.3 * (1 + z)**3 + 0.7)
    Hz = Hz / 1000.0
    return (G * M200 / (100 * Hz**2))**(1 / 3)

z_list = [0.0]

for z in z_list:
    df = pd.read_csv("random_galaxy_results_z1.txt", delim_whitespace=True)
    df = df[np.isclose(df["z"], z)]

    fgas = df["Mgas"] / df["M200"]
    fstar = df["Mstellar"] / df["M200"]
    fbary = (df["Mgas"] + df["Mstellar"]) / df["M200"]

    mass_bins = np.array([1e10, 3e10, 1e11, 3e11, 1e12, 3e12])

    success_text = []

    for low, high in zip(mass_bins[:-1], mass_bins[1:]):
        mask = (df["M200"] >= low) & (df["M200"] < high)
        n_total = mask.sum()
        n_success = df.loc[mask, "success"].sum()
        frac = 100 * n_success / n_total if n_total > 0 else np.nan

        success_text.append(
            f"{low:.0e}–{high:.0e}: {int(n_success)}/{int(n_total)} ({frac:.1f}%)"
        )

    plt.rcParams.update({
        "font.size": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
    })

    fig, axes = plt.subplots(2, 3, figsize=(7.0, 3.8))

    cmap = "viridis"
    Rmax = 3**0.5 * 3 * R_200(df["M200"], z)

    norm = LogNorm(
        vmin=(df["radius_max"] / Rmax).min(),
        vmax=(df["radius_max"] / Rmax).max()
    )

    ax = axes[0, 0]

    for y in ["Rgas", "Rstar", "Hgas", "Hstar"]:
        ax.scatter(
            df["M200"],
            df[y],
            c=df["radius_max"] / Rmax,
            cmap=cmap,
            norm=norm,
            marker="x",
            s=12,
            alpha=0.5,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"$R,H\ [{\rm kpc}]$")

    ax = axes[0, 1]

    ax.scatter(
        df["M200"],
        df["concentration"],
        c=df["radius_max"] / Rmax,
        cmap=cmap,
        norm=norm,
        marker="x",
        s=12,
        alpha=0.5,
    )

    ax.set_xscale("log")
    ax.set_ylabel(r"$c$")

    ax = axes[1, 0]

    for y in [fgas, fstar, fbary]:
        ax.scatter(
            df["M200"],
            y,
            c=df["radius_max"] / Rmax,
            cmap=cmap,
            norm=norm,
            marker="x",
            s=12,
            alpha=0.5,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{200}\ [{\rm M_\odot}]$")
    ax.set_ylabel("Mass fraction")

    ax = axes[1, 1]

    ax.scatter(
        df["M200"],
        df["Mdot"],
        c=df["radius_max"] / Rmax,
        cmap=cmap,
        norm=norm,
        marker="x",
        s=12,
        alpha=0.5,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{200}\ [{\rm M_\odot}]$")
    ax.set_ylabel(r"$\dot M_\star\ [{\rm M_\odot\,yr^{-1}}]$")

    ax = axes[1, 2]

    ax.scatter(
        df["M200"],
        df["Z0"],
        c=df["radius_max"] / Rmax,
        cmap=cmap,
        norm=norm,
        marker="x",
        s=12,
        alpha=0.5,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{200}\ [{\rm M_\odot}]$")
    ax.set_ylabel(r"$Z_0/Z_\odot$")

    axes[0, 2].axis("off")

    axes[0, 2].text(
        0.0,
        1.0,
        "Successful galaxies\n\n" + "\n".join(success_text),
        transform=axes[0, 2].transAxes,
        ha="left",
        va="top",
        fontsize=6,
    )

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cax = inset_axes(
        axes[0, 2],
        width="90%",
        height="12%",
        loc="lower center",
        borderpad=0.8,
    )

    cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label(r"$R_{\rm max}$", fontsize=6)

    plt.tight_layout()
    plt.savefig(
        "random_catalogue_scatter_z{}.png".format(z),
        dpi=300,
    )
    plt.close()