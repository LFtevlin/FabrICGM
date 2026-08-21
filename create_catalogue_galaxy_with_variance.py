import subprocess
import tempfile
import time

template_file = "param_crv-tmp.txt"
cooling_file = "./cooling/UVB_dust1_CR1_G1_shield0.hdf5"

with open(template_file) as f:
    template = f.read()

def run_validation_simulation():
    text = template
    text = text.replace("validate_galaxies=..", "validate_galaxies=true")
    text = text.replace("allow_variance=..", "allow_variance=true")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=True) as tmp:
        tmp.write(text)
        tmp.flush()

        proc = subprocess.run(
            ["./CGMake", tmp.name, cooling_file],
            capture_output=True,
            text=True,
        )

    print(proc.stdout)

    if proc.stderr:
        print(proc.stderr)

    lines = proc.stdout.strip().splitlines()
    catalogue_line = None

    for line in lines:
        if line.startswith("CATALOGUE_RESULT"):
            catalogue_line = line
            break

    if catalogue_line is None:
        return None

    values = catalogue_line.split()[1:]

    if len(values) < 14:
        return None

    return {
        "M200": values[0],
        "z": values[1],
        "success": int(values[2]),
        "radius_max": values[3],
        "concentration": values[4],
        "Z0": values[5],
        "Mgas": values[6],
        "Mstellar": values[7],
        "Rgas": values[8],
        "Rstar": values[9],
        "Hgas": values[10],
        "Hstar": values[11],
        "Mdot": values[12],
        "Rsonic": values[13],
    }

print("Searching for successful realization...")

start = time.time()
n_try = 0
good_result = None

while good_result is None:
    n_try += 1
    print(f"\nTrying realization {n_try}")

    result = run_validation_simulation()

    if result is None:
        continue

    print("RESULT:", result)

    if result["success"]:
        good_result = result
        elapsed = time.time() - start
        print("\nFound solution!")
        print(f"Realizations tried: {n_try}")
        print(f"Wall time: {elapsed:.1f} s")

text = template

def replace_parameter(text, key, value):
    lines = text.splitlines()

    for i, line in enumerate(lines):
        if line.startswith(key + "="):
            lines[i] = f"{key}={value}"

    return "\n".join(lines)

text = replace_parameter(text, "validate_galaxies", "false")
text = replace_parameter(text, "create_random_galaxy_catalogue", "false")
text = replace_parameter(text, "allow_variance", "false")

for key in [
    "M200",
    "z",
    "concentration",
    "Z0",
    "Mgas",
    "Mstellar",
    "Rgas",
    "Rstar",
    "Hgas",
    "Hstar",
    "Mdot",
    "Rsonic",
]:
    text = replace_parameter(text, key, good_result[key])

print("\nRunning final simulation with fixed parameters...")

with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=True) as tmp:
    tmp.write(text)
    tmp.flush()

    proc = subprocess.run(
        ["./CGMake", tmp.name, cooling_file],
        capture_output=True,
        text=True,
    )

print(proc.stdout)

if proc.stderr:
    print(proc.stderr)

print("Done.")