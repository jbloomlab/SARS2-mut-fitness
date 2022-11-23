library(ggplot2)
library(ggforce)

`%notin%` <- Negate(`%in%`)

# Create filename for plotting mutation data from a two-clade comparison.
input_name <- paste("results/pvalues", snakemake@wildcards[["comparisons"]], "_", snakemake@wildcards[["mutation_type"]], ".csv", sep="")

# Read in the data for plotting.
verify_analyses <- read.csv(input_name, sep=",", header=TRUE)

# Plot the corrected p-values for synonymous and non-synonymous mutations.
if("synonymous" %notin% colnames(verify_analyses)){
  verify_analyses$synonymous <- ifelse(verify_analyses$synonymous_1=="True" & verify_analyses$synonymous_2=="True", "True", "False")
}
verify_analyses$p_one <- ifelse(verify_analyses$p_value==1, "P=1", "P<1")

p <- ggplot(verify_analyses, aes(x=p_value, fill=synonymous, label=after_stat(count)))+
  geom_histogram()+
  geom_text(size=2, vjust=-0.5, stat="bin", position="stack")+
  geom_vline(data=data.frame(xint=0.05, p_one="P<1"), aes(xintercept=xint), linetype=2)+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))+
  labs(x="p-values")+
  theme_bw()

pdf(file=snakemake@output[[1]], width=12, height=10)
  
p + ggforce::facet_row(vars(p_one), scales="free", space="free")

dev.off()
