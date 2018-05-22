library(ggplot2)
library(reshape2)
library(cowplot)
setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/Silva_org/")

HMP <- read.table("HMP_Unoise_Med_silva_org.tsv", header=T, sep="\t")

HMP$Mean <- rowMeans(HMP[,c(2,3,4,5)])
HMP <- HMP[, -c(2,3,4,5)]

HMP_melt <- melt(data=HMP, measure.vars=c("Mean"), 
                 id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")

HMP_plot <- ggplot(HMP_melt, aes(x = sample, y = rel_abun, fill = Organism)) + geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("#0080ff", "#33cc33", "#ff0000")) + labs(x = "Mean Abundance Across Samples", y = "Relative Abundance") + 
  theme(axis.text.x = element_blank())
HMP_plot


Zymock <- read.table("Zymock_Unoise_Med_silva_org.tsv", header=T, sep="\t")

Zymock$Mean <- rowMeans(Zymock[, c(2,3,4)])
Zymock <- Zymock[, -c(2,3,4)]

Zymock_melt <- melt(data=Zymock, measure.vars=c("Mean"),
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
Zymock_plot <- ggplot(Zymock_melt, aes(x = sample, y = rel_abun, fill=Organism)) + geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("#3399ff", "#ff0000")) + labs(x = "Mean Abundance Across Samples", y = "Relative Abundance") + 
  theme(axis.text.x = element_blank())
Zymock_plot

grid_plot <- plot_grid(HMP_plot, Zymock_plot, labels=c("A", "B"))
grid_plot


HMP$precent <- HMP$Mean/ sum(HMP$Mean)

Zymock$precent <- Zymock$Mean / sum(Zymock$Mean)
