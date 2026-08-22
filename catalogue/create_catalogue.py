import numpy as np

def build_loglog_linear_function(x1, y1, x2, y2):
    logx1, logx2 = np.log10(x1), np.log10(x2)
    logy1, logy2 = np.log10(y1), np.log10(y2)
    k = (logy2 - logy1) / (logx2 - logx1)

    A = 10**(logy1 - k * logx1)

    def f(x):
        return A * x**k

    return f

G = 4.30091e-6  # kpc (km/s)^2 Msun^-1
H0 = 70.0       # km/s/Mpc


def R_200(M200, z):
    Hz = H0 * np.sqrt(0.3*(1+z)**3 + 0.7)  # km/s/Mpc
    Hz = Hz / 1000.0                       # km/s/kpc
    
    R = (G*M200/(100*Hz**2))**(1/3)
    return R  # kpc


## this adds a redshift dependence to the stellar mass, fits together with Behroozi+2019
def M_stars(z, M_stellar):
    nu = -0.72
    return M_stellar * (z+1)**nu

def Z_gas(z,M_star): ## in Z_solar, normalize so that ~1 in 1e10M_star at z=0, from Fig.5 Xiangcheng Ma 2015
        ## in inflow: factor 3 smaller, Fig. 11
        gamma_g = 0.35
        Z_g10 = 0.93*np.exp(-0.43*z) - 1.05

        return 10**(gamma_g * (np.log10(M_star) - 10) + Z_g10)

## both scaling relations from Diemer&Kravtsov 2013
def Rstar(M200, z):
     R200 = R_200(M200, z)
     return R200 * 0.011 

def Rgas(M200, z):
     R200 = R_200(M200, z)
     return R200 * 0.029 


## concentration parameter from Fig. 10 Dutton+2014
## z=0, x-axis: log10 M_sun/h, y-axis: log10 c
c0_Dutton = build_loglog_linear_function(10, 1.1, 12.5, 0.85) 
## z=0.5, x-axis: log10 M_sun/h, y-axis: log10 c
c05_Dutton = build_loglog_linear_function(10, 0.98, 12.5, 0.77)
## z=1, x-axis: log10 M_sun/h, y-axis: log10 c
c1_Dutton = build_loglog_linear_function(10, 0.88, 12.5, 0.69)  
## z=2, x-axis: log10 M_sun/h, y-axis: log10 c
c2_Dutton = build_loglog_linear_function(10, 0.74, 12.5, 0.585) 

## Mstellar from Behroozi et al. (2019) – UniverseMachine, Fig. 9, right
## z=0.1, x-axis M200 in M_sun, y-axis
## values on x-axis:
M200_Behroozi = np.array([2e10, 5e10, 1e11, 5e11, 1e12, 5e12])
## values on y-axis:
Mstellar_Behroozi = np.array([1e7, 8e7, 3.2e8, 8e9, 2.2e10, 8e10])

## Mstellar from Fig. 1 Dev+2024
## xray selected nearby galaxies, unit Msun
M200_Dev = [1e10, 5e10, 1e11, 5e11, 1e12, 5e12]
Mstellar_Dev = [3.2*1e7, 2.2*1e8, 1.4*1e9, 8.5*1e9, 2.9*1e10, 7*1e10]


## sSFR from Behroozi et al. (2019) – UniverseMachine, Fig. 3, right
## z=0.1, x-axis Mstellar in M_sun,
Mstellar_Behroozi_xaxis = np.array([1e8, 5e8, 1e9, 5e9, 1e10, 5e10])
## values on y-axis, sSFR in yr^-1:
Mdot_z01_Behroozi = np.array([3e-10, 2.5e-10, 2.3e-10, 1.5e-10, 1.2e-10, 4.5e-11])
## z=1
Mdot_z1_Behroozi = np.array([5e-10, 9e-10, 1e-9, 1.2e-9, 8e-10, 3e-10])
## z=2
Mdot_z2_Behroozi = np.array([1e-9, 2e-9, 2e-9, 2.5e-9, 2.5e-9, 1.5e-9])

## SFR from Spitoni 2020 Fig 1
## z=2
## stellar mass (in Msun) vs. SFR (in Msun/yr)
f2_Spitoni = build_loglog_linear_function(10**(8.5), 10**(-0.1), 10**(11.3), 10**(1.9)) 

##SFR from Aumer+2013 Fig. 4
## z=0
## stellar mass (in Msun) vs. SFR (in Msun/yr)
f0_Aumer = build_loglog_linear_function(1e9, 0.25, 3e11, 20) 

## Mbaryon = Mstellar + Mgas to M200 relation from Dev+2024, Fig. 3
## z=0
## y axis in units of M200
M_baryon_z0_Dev = np.array([1e-2, 10**(-1.75), 10**(-1.6), 10**(-1.43), 10**(-1.41), 10**(-1.58)])
## x-axis 
M200_z0_Dev = np.array([1e10, 5e10, 1e11, 5e11, 1e12, 5e12])
## M_baryon from Kulier+2019 for z=2, Fig 8
## y-axis
M_baryon_z2_Kulier = np.array([10**(8.5), 10**9, 10**(9.7), 10**(10.3), 10**(10.9), 10**(11.3)])
## x-axis
M200_z2_Kulier = np.array([1e10, 5e10, 1e11, 5e11, 1e12, 5e12])

def find_params(M200,z):

    R_gas = Rgas(M200,z)
    R_star = Rstar(M200,z)

    if M200 <= 5e10:
        H_gas = R_gas * 0.3
        H_star = R_star * 0.3
    else:
        H_gas = R_gas * 0.15
        H_star = R_star * 0.15


    log10M200 = np.log10(M200/0.67)
    z_vals = np.array([0, 0.5, 1, 2])
    c0 = c0_Dutton(log10M200)
    c05 = c05_Dutton(log10M200)
    c1 = c1_Dutton(log10M200)
    c2 = c2_Dutton(log10M200)
    c_vals = np.array([c0, c05, c1, c2])

    C = np.interp(z, z_vals, c_vals)


    M_stellar_z0 = np.interp(M200, M200_Behroozi, Mstellar_Behroozi)
    M_stellar = M_stars(z, M_stellar_z0)

    Mdot_z01 = np.interp(M_stellar, Mstellar_Behroozi_xaxis, Mdot_z01_Behroozi*Mstellar_Behroozi_xaxis)
    Mdot_z1 = np.interp(M_stellar, Mstellar_Behroozi_xaxis, Mdot_z1_Behroozi*Mstellar_Behroozi_xaxis)
    Mdot_z2 = np.interp(M_stellar, Mstellar_Behroozi_xaxis, Mdot_z2_Behroozi*Mstellar_Behroozi_xaxis)

    Mdot_array = np.array([Mdot_z01, Mdot_z1, Mdot_z2])

    z_vals = np.array([0.1, 1, 2])

    Mdot = np.interp(z, z_vals, Mdot_array)

    M_baryon_z0 = np.interp(M200, M200_z0_Dev, M200_z0_Dev*M_baryon_z0_Dev)
    M_baryon_z2 = np.interp(M200, M200_z2_Kulier, M_baryon_z2_Kulier)

    z_vals = np.array([0,2])
    M_baryon_vals = np.array([M_baryon_z0, M_baryon_z2])

    M_baryon = np.interp(z, z_vals, M_baryon_vals)
    M_gas = M_baryon - M_stellar

    Zgas = Z_gas(z,M_stellar)

    return R_gas, R_star, H_gas, H_star, 10**C, M_stellar, M_gas, Mdot, Zgas





import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Mass array
M200_array = np.logspace(np.log10(1e10), np.log10(5e12), 100)

# Redshifts
redshifts = [0, 0.5, 1, 1.5, 2]


# Store results
results = {}

for z in redshifts:

    Rgas_arr = []
    Rstar_arr = []
    Hgas_arr = []
    Hstar_arr = []
    C_arr = []
    Mgas_arr = []
    Mstar_arr = []
    Mdot_arr = []
    Zgas_arr = []

    for M200 in M200_array:

        Rgas1, Rstar1, Hgas1, Hstar1, C1, Mstar1, Mgas1, Mdot1, Zgas1 = find_params(M200, z)

        Rgas_arr.append(Rgas1)
        Rstar_arr.append(Rstar1)
        Hgas_arr.append(Hgas1)
        Hstar_arr.append(Hstar1)

        C_arr.append(C1)

        Mgas_arr.append(Mgas1)
        Mstar_arr.append(Mstar1)

        Mdot_arr.append(Mdot1)
        Zgas_arr.append(Zgas1)


    results[z] = {
        "Rgas": np.array(Rgas_arr),
        "Rstar": np.array(Rstar_arr),
        "Hgas": np.array(Hgas_arr),
        "Hstar": np.array(Hstar_arr),
        "C": np.array(C_arr),
        "Mgas": np.array(Mgas_arr),
        "Mstar": np.array(Mstar_arr),
        "Mdot": np.array(Mdot_arr),
        "Zgas": np.array(Zgas_arr),
    }



col_width = 3.5       # single column
text_width = 7.2      # double column

def set_figure_size(width=col_width, aspect=0.75):
    return (width, width * aspect)

col_width = 3.5  # inch (A&A single column)
panel_height = 1.8  # height per subplot (tune this!)

nrows = 5



plt.rcParams.update({
    "font.size": 7,
    "axes.labelsize": 7,
    "legend.fontsize": 5,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "lines.linewidth": 1.0,
})



fig, axes = plt.subplots(2, 3, figsize=(3.5*2, 2 * 1.5))

# fig, axes = plt.subplots(
#     2, 2,
#     figsize=(8, 6),
#     sharex=True
# )

colors = {
    0: "tab:blue",
    0.5: "tab:orange",
    1: "tab:green",
    1.5: "tab:purple",
    2: "tab:red"
}



# ----------------------------
# Legend handles
# ----------------------------

# redshift legend
redshift_handles = [
    Line2D([0], [0], color=colors[0], label=r"$z=0$"),
    Line2D([0], [0], color=colors[0.5], label=r"$z=0.5$"),
    Line2D([0], [0], color=colors[1], label=r"$z=1$"),
    Line2D([0], [0], color=colors[1.5], label=r"$z=1.5$"),
    Line2D([0], [0], color=colors[2], label=r"$z=2$")
]


# size line-style legend
size_handles = [
    Line2D([0], [0], color="black",
           linestyle="-", label=r"$R_{\rm gas}$"),

    Line2D([0], [0], color="black",
           linestyle="--", label=r"$R_\star$"),

    Line2D([0], [0], color="black",
           linestyle=":", label=r"$H_{\rm gas}$"),

    Line2D([0], [0], color="black",
           linestyle="-.", label=r"$H_\star$")
]


# baryon fraction line-style legend
baryon_handles = [
    Line2D([0], [0], color="black",
           linestyle="-",
           label=r"$M_{\rm gas}/M_{200}$"),

    Line2D([0], [0], color="black",
           linestyle="--",
           label=r"$M_\star/M_{200}$"),

    Line2D([0], [0], color="black",
           linestyle=":",
           label=r"$(M_{\rm gas}+M_\star)/M_{200}$")
]



# ----------------------------
# Panel 1: Rgas, Rstar, Hgas, Hstar
# ----------------------------

ax = axes[0,0]

for z in redshifts:

    ax.plot(M200_array,
            results[z]["Rgas"],
            color=colors[z])

    ax.plot(M200_array,
            results[z]["Rstar"],
            color=colors[z],
            linestyle="--")

    ax.plot(M200_array,
            results[z]["Hgas"],
            color=colors[z],
            linestyle=":")

    ax.plot(M200_array,
            results[z]["Hstar"],
            color=colors[z],
            linestyle="-.")


ax.set_xscale("log")
ax.set_yscale("log")

ax.set_ylabel(r"$R,\ H\ [{\rm kpc}]$")

ax.legend(handles=size_handles, frameon=False, loc="best")



# ----------------------------
# Panel 2: concentration
# ----------------------------

ax = axes[0,1]

for z in redshifts:

    ax.plot(M200_array,
            results[z]["C"],
            color=colors[z])


ax.set_xscale("log")

ax.set_ylabel(r"$C$")

ax.legend(handles=redshift_handles,frameon=False,
          loc="best")



# ----------------------------
# Panel 3: baryon fractions
# ----------------------------

ax = axes[1,0]


for z in redshifts:

    Mgas = results[z]["Mgas"]
    Mstar = results[z]["Mstar"]


    ax.plot(M200_array,
            Mgas/M200_array,
            color=colors[z])


    ax.plot(M200_array,
            Mstar/M200_array,
            color=colors[z],
            linestyle="--")


    ax.plot(M200_array,
            (Mgas+Mstar)/M200_array,
            color=colors[z],
            linestyle=":")



ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$M_{200}\ [M_\odot]$")
ax.set_ylabel("Mass fraction")



ax.legend(handles=baryon_handles,
          frameon=False,
          loc="best")



# ----------------------------
# Panel 4: SFR
# ----------------------------

ax = axes[1,1]


for z in redshifts:

    ax.plot(M200_array,
            results[z]["Mdot"],
            color=colors[z])


ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$M_{200}\ [M_\odot]$")
ax.set_ylabel(r"$\dot M_\star\ [{\rm M_\odot\,yr^{-1}}]$")

ax = axes[1,2]


for z in redshifts:

    ax.plot(M200_array,
            results[z]["Zgas"],
            color=colors[z])


ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"$M_{200}\ [M_\odot]$")
ax.set_ylabel(r"$Z_0/\mathrm{Z}_\odot$")


# no legend here



# ----------------------------
# Final formatting
# ----------------------------

# for ax in axes.flat:
#     ax.grid(alpha=0.3)


plt.tight_layout()
plt.savefig('./catalogue/catalogue.pdf')

import h5py


filename = "./catalogue/galaxy_catalogue.hdf5"





with h5py.File(filename, "w") as f:

    # independent axes
    f.create_dataset(
        "M200",
        data=M200_array
    )

    f.create_dataset(
        "redshift",
        data=np.array(redshifts)
    )


    # create 2D arrays (z, M200)
    Rgas = np.array([results[z]["Rgas"] for z in redshifts])
    Rstar = np.array([results[z]["Rstar"] for z in redshifts])
    Hgas = np.array([results[z]["Hgas"] for z in redshifts])
    Hstar = np.array([results[z]["Hstar"] for z in redshifts])

    C = np.array([results[z]["C"] for z in redshifts])

    Mgas = np.array([results[z]["Mgas"] for z in redshifts])
    Mstar = np.array([results[z]["Mstar"] for z in redshifts])

    Mdot = np.array([results[z]["Mdot"] for z in redshifts])
    Zgas = np.array([results[z]["Zgas"] for z in redshifts])


    # store parameters
    f.create_dataset("Rgas", data=Rgas)
    f.create_dataset("Rstar", data=Rstar)

    f.create_dataset("Hgas", data=Hgas)
    f.create_dataset("Hstar", data=Hstar)

    f.create_dataset("Concentration", data=C)

    f.create_dataset("Mgas", data=Mgas)
    f.create_dataset("Mstar", data=Mstar)

    f.create_dataset("Mdot", data=Mdot)
    f.create_dataset("Zgas", data=Zgas)


print(f"Saved {filename}")