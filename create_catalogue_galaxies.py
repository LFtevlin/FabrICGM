import numpy as np
from tempfile import NamedTemporaryFile
import subprocess
import h5py
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import glob
import shutil

z_grid = [0]
n_random = 5000
N_workers = 16

with open("./param_cc-tmp.txt") as f:
    template = f.read()

cooling_file = "./cooling/UVB_dust1_CR1_G1_shield0.hdf5"

def run_simulation(job):
    M200, z, realization = job

    text = template
    text = text.replace("M200=..", f"M200={M200}")
    text = text.replace("z=..", f"z={z}")
    text = text.replace("NAME", f"M{M200:.3e}_z{z:.2f}_r{realization}")

    with NamedTemporaryFile(mode="w", suffix=".txt", delete=True) as tmp:
        tmp.write(text)
        tmp.flush()

        try:
            proc = subprocess.run(
                ["./FabriCGM", tmp.name, cooling_file],
                capture_output=True,
                text=True,
                timeout=150
            )
        except subprocess.TimeoutExpired:
            return {
                "M200": M200,
                "z": z,
                "success": 0,
                "radius_max": 0,
            }

    for line in proc.stdout.splitlines():
        if line.startswith("CATALOGUE_RESULT"):
            values = line.split()

            return {
                "M200": float(values[1]),
                "z": float(values[2]),
                "success": int(values[3]),
                "radius_max": float(values[4]),
                "concentration": float(values[5]),
                "Z0": float(values[6]),
                "Mgas": float(values[7]),
                "Mstellar": float(values[8]),
                "Rgas": float(values[9]),
                "Rstar": float(values[10]),
                "Hgas": float(values[11]),
                "Hstar": float(values[12]),
                "Mdot": float(values[13]),
                "Rsonic": float(values[14]),
                "Mcgm": float(values[15]),
            }

    return {
        "M200": M200,
        "z": z,
        "success": 0,
        "radius_max": 0,
    }

jobs = []

for z in z_grid:
    for i in range(n_random):
        M200 = 10**np.random.uniform(np.log10(1e10), np.log10(1e12))
        jobs.append((M200, z, i))

print(f"Total simulations: {len(jobs)}")

results = []
output_file = "./catalogue/random_galaxy_results2.txt"

with open(output_file, "w") as f:
    f.write(
        "M200 z success radius_max concentration Z0 "
        "Mgas Mstellar Rgas Rstar Hgas Hstar Mdot Rsonic Mcgm\n"
    )

with ProcessPoolExecutor(max_workers=N_workers) as executor:
    futures = [executor.submit(run_simulation, job) for job in jobs]

    with open(output_file, "a") as outfile:
        for i, future in enumerate(as_completed(futures)):
            result = future.result()

            print(
                f"[{i+1}/{len(jobs)}] "
                f"M200={result['M200']:.2e}, "
                f"z={result['z']:.2f}, "
                f"success={bool(result['success'])}, "
                f"Rmax={result['radius_max']:.2f}"
            )

            outfile.write(
                f"{result['M200']} "
                f"{result['z']} "
                f"{result['success']} "
                f"{result['radius_max']} "
                f"{result.get('concentration', np.nan)} "
                f"{result.get('Z0', np.nan)} "
                f"{result.get('Mgas', np.nan)} "
                f"{result.get('Mstellar', np.nan)} "
                f"{result.get('Rgas', np.nan)} "
                f"{result.get('Rstar', np.nan)} "
                f"{result.get('Hgas', np.nan)} "
                f"{result.get('Hstar', np.nan)} "
                f"{result.get('Mdot', np.nan)} "
                f"{result.get('Rsonic', np.nan)} "
                f"{result.get('Mcgm', np.nan)}\n"
            )

            if (i + 1) % 100 == 0:
                for filename in glob.glob("out_M*.txt"):
                    try:
                        os.remove(filename)
                    except FileNotFoundError:
                        pass