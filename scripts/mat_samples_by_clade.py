import os

import pandas as pd


os.makedirs(snakemake.output.subdir, exist_ok=True)

for clade, df in pd.read_csv(snakemake.input.csv).groupby("nextstrain_clade"):
    clade = clade.split()[0]
    outfile = os.path.join(snakemake.output.subdir, f"{clade}.csv")
    print(f"Writing {len(df)} samples for {clade} to {outfile}")
    df.to_csv(outfile, index=False)
