import pandas as pd


(
    pd.read_csv(
        snakemake.input.vcf,
        sep="\t",
        comment="#",
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
    )
    .rename(columns={"POS": "site"})
    [["site"]]
    .to_csv(snakemake.output.csv)
)
