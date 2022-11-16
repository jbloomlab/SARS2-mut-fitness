library(dplyr)

# Read data using access to the Fred Hutch server.
SARS_mutation <- read.csv(snakemake@input[[1]], sep=",", header=TRUE)

# Include data from certain countries AND exclude branches with very high mutation rates and questionable mutations.
filtered_SARS_mutation <- subset(SARS_mutation, subset==snakemake@params@subset & exclude=="False" & expected_count>=snakemake@params@min_expected)

# Make a list of the unique clades in the dataset.
clades <- unique(filtered_SARS_mutation$clade)

# Filter columns to make workflow easier to interpret.
if(snakemake@wildcards[["nt"]]=="nt"){
   nt_group <- filtered_SARS_mutation[,c("clade","nt_mutation","clade_founder_nt","nt_site","synonymous","four_fold_degenerate","gene","actual_count","expected_count")]
}

# Group by amino acid mutation and add counts for identical amino acid mutations.
if(snakemake@wildcards[["aa"]]=="aa"){
   aa_group <- as.data.frame(filtered_SARS_mutation %>% group_by(clade, aa_mutation, clade_founder_aa, codon_site, mutant_aa, synonymous, gene) %>% summarise(actual_count=sum(actual_count), expected_count=sum(expected_count)))
}

#
filtered_SARS_mutation <-  








# Calculate the ratio of actual counts (i.e., observed) to expected counts (i.e., expected) with a pseudocount.
filtered_SARS_mutation$observed_expected <- (filtered_SARS_mutation$actual_count + snakemake@input[[2]])/(filtered_SARS_mutation$expected_count + snakemake@input[[2]])

# Round expected counts because contingency table analyses (e.g., Fisher's Exact test) can only analyze integers.
filtered_SARS_mutation$rounded_expected <- round(filtered_SARS_mutation$expected_count)

# Create two data frames with data for clade 1 and clade 2.
clade1 <- filtered_SARS_mutation[which(filtered_SARS_mutation$clade==clades[1]),]
clade2 <- filtered_SARS_mutation[which(filtered_SARS_mutation$clade==clades[2]),]

# Append '1' or '2' to the column names so that the two data frames can be easily combined.
colnames(clade1)[colnames(clade1) %notin% merge_column] <- paste(colnames(clade1)[colnames(clade1) %notin% merge_column], "1", sep="_")
colnames(clade2)[colnames(clade2) %notin% merge_column] <- paste(colnames(clade2)[colnames(clade2) %notin% merge_column], "2", sep="_")

# Merge the two data frames by mutation so that each row represents one comparison.
comparisons <- merge(clade1, clade2, by=merge_column)

# Determine if any of the columns from 'clade1' and 'clade2' are the same.
col_names <- colnames(filtered_SARS_mutation)[colnames(filtered_SARS_mutation) %notin% merge_column]
is_equal <- data.frame(sapply(col_names, FUN=function(c){identical(comparisons[,paste(c[1], 1, sep="_")], comparisons[,paste(c[1], 2, sep="_")])}))
col_equal <- row.names(is_equal)[which(is_equal[1]=="TRUE")]

# Rename one of the duplicate columns and remove the other.
colnames(comparisons)[colnames(comparisons) %in% paste(col_equal, "1", sep="_")] <- col_equal
comparisons[,paste(col_equal, "2", sep="_")] <- NULL

# Iterate through rows and compare clade 1 observed and expected values to clade 2 observed and expected values.
p_values <- data.frame(p_value=apply(comparisons, 1, FUN=function(c) {fisher.test(rbind(c(as.numeric(c["actual_count_1"]),as.numeric(c["rounded_expected_1"])),c(as.numeric(c["actual_count_2"]),as.numeric(c["rounded_expected_2"]))))$p.value}))

# Combine the p-value list with the comparison list.
comparison_p_values <- cbind(comparisons, p_values)

# Using a false discovery rate correction, adjust the p-values by the number of hypothesis tests.
comparison_p_values$p_value_corrected <- p.adjust(comparison_p_values$p_value, method="fdr", n=nrow(comparison_p_values))

# Calculate fold-change for plotting volcano plots.
comparison_p_values$fold_change <- log2(comparison_p_values$observed_expected_1/comparison_p_values$observed_expected_2)
