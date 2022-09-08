import os

import pandas as pd


min_clade_samples = snakemake.params.min_clade_samples

os.makedirs(snakemake.output.subdir, exist_ok=True)

for clade, df in pd.read_csv(snakemake.input.csv).groupby("nextstrain_clade"):
    clade = clade.split()[0]
    if len(df) < min_clade_samples:
        print(f"Ignoring {clade=} as it has only {len(df)} samples")
    else:
        outfile = os.path.join(snakemake.output.subdir, f"{clade}.txt")
        print(f"Writing {len(df)} samples for {clade} to {outfile}")
        df["sample"].to_csv(outfile, index=False, header=False)
