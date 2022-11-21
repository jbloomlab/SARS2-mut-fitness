library(ggplot2)
library(RColorBrewer)

`%notin%` <- Negate(`%in%`)

# Specify column for plotting mutation data (nt_mutation or aa_mutation).
merge_column <- paste(snakemake@wildcards[["mutation_type"]], "mutation", sep="_")

# Create filename for plotting mutation data from a two-clade comparison.
input_name <- paste("results/pvalues", snakemake@wildcards[["comparisons"]], "_", snakemake@wildcards[["mutation_type"]], ".csv", sep="")

# Read in the data for plotting.
comparison_p_values <- read.csv(input_name, sep=",", header=TRUE)

## After dropping P = 1, plot the corrected p-values. ##
volcano <- comparison_p_values[which(comparison_p_values$p_value_corrected!=1),]
if("gene" %notin% colnames(volcano)){
  volcano$gene <- ifelse(volcano$gene_1==volcano$gene_2, volcano$gene_1, paste(volcano$gene_1, volcano$gene_2, sep="_"))
}

if(merge_column == "nt_mutation"){
	color <- colorRampPalette(brewer.pal(8,"Set1"))(length(unique(volcano$gene)))

	plot <- ggplot(volcano, aes(x=fold_change, y=-log10(p_value_corrected), label=get(merge_column), color=gene))+
  	  geom_point()+
  	  geom_text(color="black", size=2, nudge_y=0.2)+
  	  geom_hline(yintercept=-log10(0.05), col="red", linetype=2)+
  	  scale_color_manual(values=color, name="mutation")+
  	  theme_bw()+
  	  theme(legend.position="bottom")
}

## Aggregate all mutations that do not include the spike protein to 'non-spike'.
if(merge_column == "aa_mutation"){
	spike <- volcano$gene[grep("S", volcano$gene)]
	volcano$spike <- ifelse(volcano$gene %in% spike, volcano$gene, "non-spike")
	non_spike <- which(levels(as.factor(volcano$spike))=="non-spike")
	color <- colorRampPalette(brewer.pal(8,"Set1"))(length(unique(volcano$spike)))
	color[non_spike] <- "grey"

	plot <- ggplot(volcano, aes(x=fold_change, y=-log10(p_value_corrected), label=get(merge_column), color=spike))+
	  geom_point()+
	  geom_text(color="black", size=2, nudge_y=0.2)+
	  geom_hline(yintercept=-log10(0.05), col="red", linetype=2)+
	  scale_color_manual(values=color, name="mutation")+
	  theme_bw()
}

pdf(file=snakemake@output[[1]])

plot

dev.off()
