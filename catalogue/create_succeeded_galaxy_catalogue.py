import pandas as pd


# ----------------------------------------
# Input / output files
# ----------------------------------------

input_file = "random_galaxy_results_z2.txt"
output_file = "successful_galaxies_sorted_z2.txt"


# ----------------------------------------
# Read catalogue
# ----------------------------------------

df = pd.read_csv(
    input_file,
    delim_whitespace=True
)


# ----------------------------------------
# Select successful galaxies
# ----------------------------------------

df_success = df[df["success"] == 1].copy()


# ----------------------------------------
# Sort by M200
# ----------------------------------------

df_success = df_success.sort_values(
    by="M200"
)


# ----------------------------------------
# Write new catalogue
# ----------------------------------------

df_success.to_csv(
    output_file,
    sep=" ",
    index=False,
    float_format="%.6e"
)


print(
    f"Found {len(df_success)} successful galaxies"
)

print(
    f"Written to {output_file}"
)