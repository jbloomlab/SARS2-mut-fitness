library(ggplot2)
library(RColorBrewer)
library(grid)

`%notin%` <- Negate(`%in%`)

# Specify column for plotting mutation data (nt_mutation or aa_mutation).
merge_column <- paste(snakemake@wildcards[["mutation_type"]], "mutation", sep="_")

# Create filename for plotting mutation data from a two-clade comparison.
input_name <- paste("results/pvalues", snakemake@wildcards[["comparisons"]], "_", snakemake@wildcards[["mutation_type"]], ".csv", sep="")

# Read in the data for plotting.
comparison_p_values <- read.csv(input_name, sep=",", header=TRUE)

## After dropping P = 1, plot the corrected p-values. ##
volcano <- comparison_p_values[which(comparison_p_values$p_value_corrected!=1),]
color <- colorRampPalette(brewer.pal(8,"Set1"))(length(unique(volcano$gene)))

neg_clade <- textGrob(unique(volcano$clade_2), gp=gpar(fontsize=13, fontface="bold"))
pos_clade <- textGrob(unique(volcano$clade_1), gp=gpar(fontsize=13, fontface="bold"))

na <- ggplot(volcano, aes(x=fold_change, y=-log10(p_value_corrected), label=get(merge_column), color=gene))+
  geom_point()+
  geom_text(color="black", size=2, nudge_y=0.2)+
  geom_hline(yintercept=-log10(0.05), col="red", linetype=2)+
  scale_color_manual(values=color, name="Mutation")+
  theme_bw()+
  labs(x="Log2 Fold Change", y="-log10 FDR Corrected P-Value")

x <- layer_scales(na)$x$range$range
y <- layer_scales(na)$y$range$range

pdf(file=snakemake@output[[1]])

if(is.null(x)){
  na
} else {
  p <- na + coord_cartesian(xlim=c(-max(abs(x)), max(abs(x))))
  x_mod <- p$coordinates$limits[[1]]
  
  p + annotation_custom(neg_clade,xmin=x_mod[1]*0.7,xmax=x_mod[1]*0.7,ymin=y[2]*0.95,ymax=y[2]*0.95)+ 
      annotation_custom(pos_clade,xmin=x_mod[2]*0.7,xmax=x_mod[2]*0.7,ymin=y[2]*0.95,ymax=y[2]*0.95)
}

dev.off()
