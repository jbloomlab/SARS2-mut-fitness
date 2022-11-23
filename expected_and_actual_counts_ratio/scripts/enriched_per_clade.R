library(ggplot2)
library(RColorBrewer)

# Create the column to label points with.
mutation <- paste(snakemake@wildcards[["mutation_type"]], "mutation", sep="_")

# Create filename for reading in all mutation data.
input_name <- paste("results/pvalues", "_", snakemake@wildcards[["mutation_type"]], ".csv", sep="")

# Read in the data for plotting.
combined_mutations <- read.csv(input_name, sep=",", header=TRUE)

# Remove corrected p-values equal to 1.
filtered_mutations <- combined_mutations[which(combined_mutations$p_value_corrected != 1),]

# Filter for mutations enriched in each clade (positive fold-change is enriched in clade 1 and negative fold-change is enriched in clade 2).
clade1 <- filtered_mutations[which(filtered_mutations$fold_change > 0),c("clade_1","p_value_corrected","fold_change","gene", mutation)]
clade2 <- filtered_mutations[which(filtered_mutations$fold_change < 0),c("clade_2","p_value_corrected","fold_change","gene", mutation)]

# Change the clade column name so the two data frames have the same column names.
colnames(clade1)[1] <- "clade"
colnames(clade2)[1] <- "clade"

# Combine the clade 1 and clade 2 data frames for plotting.
combined_data <- rbind(clade1, clade2)

# Assign colors to each possible gene in the data set.
all_colors <- data.frame(color=colorRampPalette(brewer.pal(8,"Set1"))(length(unique(filtered_mutations$gene))), gene=levels(as.factor(filtered_mutations$gene)))

# Only plot data for one clade.
plot_data <- combined_data[which(combined_data$clade == snakemake@wildcards[["clades"]]),]

# Filter the assigned colors for the genes in the plotting data set.
color <- all_colors[all_colors$gene %in% unique(plot_data$gene), "color"]

# Generate plot showing the enriched mutations for each clade using all relevant comparisons.
p <- ggplot(plot_data, aes(x=abs(fold_change), y=-log10(p_value_corrected), color=gene, label=get(mutation)))+
  geom_point()+
  scale_color_manual(values=color, name="Mutation")+
  labs(x="Log2 Fold-Change", y="-log10 FDR Corrected P-Value")+
  theme_bw()

y <- layer_scales(p)$y$range$range

pdf(file=snakemake@output[[1]], width=12, height=10)

p + geom_text(color="black", size=3, nudge_y=diff(y)/50) 

dev.off()
