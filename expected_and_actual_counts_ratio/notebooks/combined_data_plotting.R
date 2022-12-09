library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(reshape)

data <- read.csv("/home/kfeldman/kfeldmann/2022/SARS2-mut-rates/expected_and_actual_counts_ratio/results/pvalues_aa.csv")

# Proportion of significantly enriched mutations for all clade comparisons
data$significant <- ifelse(data$p_value_corrected<=0.05, 1, 0)
temp <- data %>% group_by(clade_1, clade_2) %>% summarise(Proportion=sum(significant)/length(significant))
ggplot(temp, aes(x=clade_1, y=clade_2, fill=Proportion))+
  geom_tile()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  labs(x="Clade 1", y="Clade 2")

# Proportion of significantly enriched mutations for ORF8 and S across clades
data$significant <- ifelse(data$p_value_corrected<=0.05, 1, 0)
data$enriched_clade <- ifelse(data$fold_change>0, data$clade_1, data$clade_2)
temp <- as.data.frame(data %>% group_by(enriched_clade, gene) %>% summarise(Proportion=sum(significant)/length(significant)))

all_colors <- data.frame(color=colorRampPalette(brewer.pal(8,"Set1"))(length(unique(temp$gene))), gene=levels(as.factor(temp$gene)))
all_colors[,"color"] <- "black"
all_colors[which(all_colors$gene=="S"),"color"] <- "orange"
all_colors[which(all_colors$gene=="ORF8"),"color"] <- "skyblue"
ggplot(temp, aes(x=enriched_clade, y=Proportion, color=gene))+
  geom_point(size=4)+
  scale_color_manual(values=all_colors$color, name="Gene")+
  theme_bw()+
  labs(x="Clade")

# Proportion of significantly enriched mutations for each gene (clades combined; highlighting spike)
data$significant <- ifelse(data$p_value_corrected<=0.05, 1, 0)
temp <- as.data.frame(data %>% group_by(gene) %>% summarise(Proportion=sum(significant)/length(significant)))
temp <- temp[which(temp$Proportion>0),]
temp$gene <- factor(temp$gene, levels = temp[rev(order(temp$Proportion)),"gene"])

colors <- rep("navyblue", length(temp$gene))
colors[2] <- "red"
p <- ggplot(temp, aes(x=gene, y=Proportion, color=gene))+
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values = colors)+
  theme(axis.text=element_text(size=8, angle=45, vjust=0.5, hjust=1))+
  labs(x="Gene")

p + theme(legend.position="none")

# Minimum q-values for synonymous mutations by clade: mutations near deletions are colored and labeled (only spike protein)
data$enriched_clade <- ifelse(data$fold_change>0, data$clade_1, data$clade_2)
filtered <- data[which(data$gene=="S" & data$synonymous=="True"),c("enriched_clade","aa_mutation","p_value_corrected")]
temp <- as.data.frame(filtered %>% group_by(enriched_clade, aa_mutation) %>% summarise(average_p_value=min(p_value_corrected)))
temp$label <- ifelse(temp$aa_mutation=="P25P" | temp$aa_mutation=="V70V" | temp$aa_mutation=="Y145Y", temp$aa_mutation, "")

all_colors <- data.frame(color=colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(temp$aa_mutation))), aa_mutation=levels(as.factor(temp$aa_mutation)))
all_colors[,"color"] <- "black"
all_colors[which(all_colors$aa_mutation=="P25P"),"color"] <- "limegreen"
all_colors[which(all_colors$aa_mutation=="V70V"),"color"] <- "dodgerblue"
all_colors[which(all_colors$aa_mutation=="Y145Y"),"color"] <- "magenta"
p <- ggplot(temp, aes(x=enriched_clade, y=-log10(average_p_value), color=aa_mutation, label=label))+
  geom_jitter(position = position_jitter(seed = 1))+
  geom_text(size=2, position = position_jitter(seed = 1), hjust=0, vjust=0)+
  scale_color_manual(values=all_colors$color, name="Synonymous")+
  theme_bw()+
  labs(x="Clade", y="-log10(minimum q-value)")+
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2)

p + theme(legend.position="none")

# Significantly enriched mutations by clade: only T95I labeled and colored by comparison clade (only spike protein)
data$enriched_clade <- ifelse(data$fold_change>0, data$clade_1, data$clade_2)
data$other_clade <- ifelse(data$fold_change<0, data$clade_1, data$clade_2)
filtered <- data[which(data$gene=="S"),c("enriched_clade","other_clade","aa_mutation","p_value_corrected","synonymous")]
filtered$label <- ifelse(filtered$aa_mutation=="T95I", filtered$aa_mutation, "")
filtered$other_clade <- ifelse(filtered$aa_mutation=="T95I", filtered$other_clade, "")

all_colors <- data.frame(color=colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(filtered$other_clade))), other_clade=levels(as.factor(filtered$other_clade)))
all_colors[which(all_colors$other_clade==""),"color"] <- "black"
all_colors[which(all_colors$other_clade=="21J"),"color"] <- "dodgerblue"
p <- ggplot(filtered, aes(x=enriched_clade, y=-log10(p_value_corrected), color=other_clade, label=label))+
  geom_jitter(position = position_jitter(seed = 1, width = 0.2))+
  geom_text(size=2, position = position_jitter(seed = 1, width = 0.2), hjust=0, vjust=0)+
  theme_bw()+
  labs(x="Clade", y="-log10(q-value)")+
  scale_color_manual(values=all_colors$color)+
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2)

p

# 10 minimum q-values for each clade: color ordered by location along the genome (only spike protein)
data$enriched_clade <- ifelse(data$fold_change>0, data$clade_1, data$clade_2)
filtered <- data[which(data$gene=="S"),c("enriched_clade","aa_mutation","p_value_corrected")]
temp <- as.data.frame(filtered %>% group_by(enriched_clade, aa_mutation) %>% summarise(minimum_p_value=min(p_value_corrected)))
clades <- unique(temp$enriched_clade)
df <- data.frame()
for(c in clades){
  x <- temp[which(temp$enriched_clade==c),]
  y <- x[order(x$minimum_p_value),]
  z <- cbind(head(y, n=10),rep(1:1))
  df <- rbind(df, z)
}
df$site <- as.numeric(str_sub(df$aa_mutation,2,-2))
site_ordered <- df[order(df$site),]
df$aa_mutation <- factor(df$aa_mutation, levels=unique(site_ordered$aa_mutation))
all_colors <- data.frame(color=colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(df$aa_mutation))), aa_mutation=levels(df$aa_mutation))
p <- ggplot(df[order(df$site),], aes(x=enriched_clade, y=rep(1:1), fill=aa_mutation))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=all_colors$color, name="Mutation")+
  theme_bw()+
  labs(x="Clade", y="")
p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

# Minimum q-values for all mutations within the receptor binding domain (only spike protein)
data$enriched_clade <- ifelse(data$fold_change>0, data$clade_1, data$clade_2)
subset <- data[which(data$gene=="S"),c("enriched_clade","aa_mutation","p_value_corrected")]
subset$site <- as.numeric(str_sub(subset$aa_mutation,2,-2))
filtered <- subset[which(subset$site>=331 & subset$site<=530),]
temp <- as.data.frame(filtered %>% group_by(enriched_clade, aa_mutation, site) %>% summarise(minimum_p_value=min(p_value_corrected)))
temp$label <- ifelse(temp$minimum_p_value<=0.05, temp$aa_mutation, "")
site_ordered <- temp[order(temp$site),]
temp$aa_mutation <- factor(temp$aa_mutation, levels=unique(site_ordered$aa_mutation))

all_colors <- data.frame(color=colorRampPalette(brewer.pal(12,"Paired"))(length(unique(temp$aa_mutation))), aa_mutation=levels(as.factor(temp$aa_mutation)))
p <- ggplot(temp, aes(x=enriched_clade, y=-log10(minimum_p_value), color=aa_mutation, label=label))+
  geom_jitter(position = position_jitter(seed = 1))+
  geom_text(size=2, position = position_jitter(seed = 1), hjust=0, vjust=0)+
  scale_color_manual(values=all_colors$color, name="Mutation")+
  theme_bw()+
  labs(x="Clade", y="-log10(minimum q-value)")+
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2)

p + theme(legend.position="none")

# # Minimum q-values for all mutations within the receptor binding domain: color RBD sites with strong antibody escape (only spike protein)
num <- c(346,378,383,384,417,444,445,446,447,448,449,450,452,455,456,472,473,475,483,484,485,486,487,489,490,493,494,496,498,500,504)

data$enriched_clade <- ifelse(data$fold_change>0, data$clade_1, data$clade_2)
subset <- data[which(data$gene=="S"),c("enriched_clade","aa_mutation","p_value_corrected")]
subset$site <- as.numeric(str_sub(subset$aa_mutation,2,-2))
filtered <- subset[which(subset$site>=331 & subset$site<=530),]
temp <- as.data.frame(filtered %>% group_by(enriched_clade, aa_mutation, site) %>% summarise(minimum_p_value=min(p_value_corrected)))
temp$label <- ifelse(temp$minimum_p_value<=0.05, temp$aa_mutation, "")
site_ordered <- temp[order(temp$site),]
temp$aa_mutation <- factor(temp$aa_mutation, levels=unique(site_ordered$aa_mutation))

all_colors <- data.frame(color=colorRampPalette(brewer.pal(12,"Paired"))(length(unique(temp$aa_mutation))), aa_mutation=levels(as.factor(temp$aa_mutation)))
escape <- temp[which(temp$site %in% num),"aa_mutation"]
all_colors[which(all_colors$aa_mutation %in% escape),"color"] <- "black"

p <- ggplot(temp, aes(x=enriched_clade, y=-log10(minimum_p_value), color=aa_mutation, label=label))+
  geom_jitter(position = position_jitter(seed = 1))+
  geom_text(size=2, position = position_jitter(seed = 1), hjust=0, vjust=0)+
  scale_color_manual(values=all_colors$color, name="Mutation")+
  theme_bw()+
  labs(x="Clade", y="-log10(minimum q-value)")+
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2)

p + theme(legend.position="none")

# Generate random proportions of significant synonymous mutations and compare to real proportion
data$synonymous_num <- ifelse(data$synonymous=="True", 1, 0)
sig <- data[which(data$p_value_corrected<=0.05),]
real_sig <- sum(sig$synonymous_num)/nrow(sig)

samples <- vector()
reps <- rep(1:10000)
for(r in reps){
  temp <- data.frame(p_value_corrected=data$p_value_corrected, synonymous_num=sample(data$synonymous_num))
  sig <- temp[which(temp$p_value_corrected<=0.05),]
  sample_sig <- sum(sig$synonymous_num)/nrow(sig)
  samples[r] <- sample_sig
}
ggplot()+
  aes(samples) + geom_histogram()+ 
  geom_vline(xintercept = real_sig, color = "red")+
  theme_bw()+
  labs(x="proportion")

p_value <- sum(samples <= real_sig)/length(samples)

# Observed and expected counts for the T95I mutation in the spike protein
SARS_mutation <- read.csv("/home/kfeldman/kfeldmann/2022/SARS2-mut-rates/results/expected_vs_actual_mut_counts/expected_vs_actual_mut_counts.csv", sep=",", header=TRUE)
country <- "all"
min_expected <- 5
filtered_SARS_mutation <- subset(SARS_mutation, subset==country & exclude=="False" & expected_count>=min_expected)
t95i <- filtered_SARS_mutation[which(filtered_SARS_mutation$aa_mutation=="T95I"),]
colnames(t95i)[6] <- "Expected"
colnames(t95i)[7] <- "Observed"
plot_data <- melt(t95i[,c("clade","Observed","Expected")], id=c("clade"))
ggplot(plot_data, aes(x=clade, y=value, fill=variable))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_bw()+
  labs(x="Clade", y="Count")
