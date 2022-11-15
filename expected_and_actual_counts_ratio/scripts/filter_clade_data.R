library(dplyr)

# Read data using access to the Fred Hutch server.
SARS_mutation <- read.csv(snakemake@input[[1]], sep=",", header=TRUE)

# Include data from certain countries AND exclude branches with very high mutation rates and questionable mutations.
filtered_SARS_mutation <- subset(SARS_mutation, subset==snakemake@input[[2]] & exclude=="False" & expected_count>=snakemake@input[[3]])

# Make a list of the unique clades in the dataset.
clades <- unique(filtered_SARS_mutation$clade)

## NUCLEOTIDE MUTATION
# Filter columns to make workflow easier to interpret.
nt_group <- filtered_SARS_mutation[,c("clade","nt_mutation","clade_founder_nt","nt_site","synonymous","four_fold_degenerate","gene","actual_count","expected_count")]

# Save individual CSV files for each clade with the filtered data.
sapply(clades, FUN=function(c){ write.csv(nt_group[which(nt_group$clade==c[1]),], file=paste("results/", c[1], "_nt.csv", sep="")) })

## AMINO ACID MUTATION
# Group by amino acid mutation and add counts for identical amino acid mutations.
aa_group <- as.data.frame(filtered_SARS_mutation %>% group_by(clade, aa_mutation, clade_founder_aa, codon_site, mutant_aa, synonymous, gene) %>% summarise(actual_count=sum(actual_count), expected_count=sum(expected_count)))

# Save individual CSV files for each clade with the filtered data.
sapply(clades, FUN=function(c){ write.csv(aa_group[which(aa_group$clade==c[1]),], file=paste("results/", c[1], "_aa.csv", sep="")) })
