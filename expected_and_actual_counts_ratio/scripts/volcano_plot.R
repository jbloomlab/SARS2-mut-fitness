library(stringr)
library(ggplot2)
library(RColorBrewer)
library(grid)

`%notin%` <- Negate(`%in%`)

# Determine what the clades being compared are.
clades <- str_split(snakemake@wildcards[["comparisons"]], "_", 2)[[1]]

# Specify column for plotting mutation data (nt_mutation or aa_mutation).
merge_column <- paste(snakemake@wildcards[["mutation_type"]], "mutation", sep="_")

# Create filename for plotting mutation data from a two-clade comparison.
input_name <- paste("results/pvalues", snakemake@wildcards[["comparisons"]], "_", snakemake@wildcards[["mutation_type"]], ".csv", sep="")

# Read in the data for plotting.
comparison_p_values <- read.csv(input_name, sep=",", header=TRUE)

## After dropping P = 1, plot the corrected p-values. ##

all_colors <- data.frame(color=colorRampPalette(brewer.pal(8,"Set1"))(length(unique(comparison_p_values$gene))), gene=levels(as.factor(comparison_p_values$gene)))
volcano <- comparison_p_values[which(comparison_p_values$p_value_corrected!=1),]
color <- all_colors[all_colors$gene %in% unique(volcano$gene), "color"]

neg_clade <- textGrob(clades[2], gp=gpar(fontsize=20, fontface="bold"))
pos_clade <- textGrob(clades[1], gp=gpar(fontsize=20, fontface="bold"))

na <- ggplot(volcano, aes(x=fold_change, y=-log10(p_value_corrected), label=get(merge_column), color=gene))+
  geom_point()+
  geom_hline(yintercept=-log10(0.05), col="red", linetype=2)+
  scale_color_manual(values=color, name="Mutation")+
  theme_bw()+
  labs(x="Log2 Fold-Change", y="-log10 FDR Corrected P-Value")

x <- layer_scales(na)$x$range$range
y <- layer_scales(na)$y$range$range

pdf(file=snakemake@output[[1]], width=12, height=10)

if(is.null(x)){
  na
} else {
  p <- na + coord_cartesian(xlim=c(-max(abs(x)), max(abs(x)))) + geom_text(color="black", size=3, nudge_y=diff(y)/50)
  x_mod <- p$coordinates$limits[[1]]

  p + annotation_custom(neg_clade,xmin=x_mod[1]*0.7,xmax=x_mod[1]*0.7,ymin=y[2]*0.95,ymax=y[2]*0.95)+
      annotation_custom(pos_clade,xmin=x_mod[2]*0.7,xmax=x_mod[2]*0.7,ymin=y[2]*0.95,ymax=y[2]*0.95)
}

dev.off()
