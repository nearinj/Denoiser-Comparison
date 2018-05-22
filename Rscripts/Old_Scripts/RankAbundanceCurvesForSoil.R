#ASV relative abundance curve

setwd("~jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/")

library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)

Dada2_table <- read.table("dada2_fixed.tsv", sep="\t", header=F)




Deblur_table <- read.table("deblur_fixed.tsv", sep="\t", header=F)


Unoise_table <- read.table("unoise_fixed.tsv", sep="\t", header=F)


Open_table <- read.table("open-ref/open_fixed.tsv", sep="\t", header=F)


Dada2_table$Mean <- rowMeans(Dada2_table[ , -c(1, 66)])
Deblur_table$Mean <- rowMeans(Deblur_table[,-c(1, 66)])
Unoise_table$Mean <- rowMeans(Unoise_table[,-c(1, 66)])
Open_table$Mean <- rowMeans(Open_table[,-c(1, 66)])



Dada2Mean <- Dada2_table[,c(1, 67)]
Dada2Mean$Mean <- Dada2Mean$Mean/sum(Dada2Mean$Mean)
DeblurMean <- Deblur_table[,c(1, 67)]
DeblurMean$Mean <- DeblurMean$Mean/sum(DeblurMean$Mean)
UnoiseMean <- Unoise_table[,c(1,67)]
UnoiseMean$Mean <- UnoiseMean$Mean/sum(UnoiseMean$Mean)
OpenMean <- Open_table[,c(1,67)]
OpenMean$Mean <- OpenMean$Mean/sum(OpenMean$Mean)

Dada2test <- arrange(Dada2Mean, desc(Mean))
DeblurMean <- arrange(DeblurMean, desc(Mean))
UnoiseMean <- arrange(UnoiseMean, desc(Mean))
OpenMean <- arrange(OpenMean, desc(Mean))

plot(Dada2test$Mean)
plot(DeblurMean$Mean)
plot(UnoiseMean$Mean)
plot(OpenMean$Mean)

Dada2test$group <- "DADA2"
DeblurMean$group <- "Deblur"
UnoiseMean$group <- "UNOISE3"
OpenMean$group <- "Open"
Dada2test$order <- seq.int(nrow(Dada2test))
DeblurMean$order <- seq.int(nrow(DeblurMean))
UnoiseMean$order <- seq.int(nrow(UnoiseMean))
OpenMean$order <- seq.int(nrow(OpenMean))

dat = rbind(Dada2test, DeblurMean, UnoiseMean, OpenMean)

ASVs <- ggplot(data= dat, aes(x = order, y = log(Mean), color = group)) + geom_line() + labs(x="ASV Rank", y="log(Relative Abundance)") + theme_classic()
ASVs

TaxaDada2 <- Dada2_table[, c(66, 67)]
colnames(TaxaDada2)[1] <- "Taxa"

colTaxaDada2 <- ddply(TaxaDada2, .(Taxa), summarize, Sum=sum(Mean))
colTaxaDada2$Sum <- colTaxaDada2$Sum/ sum(colTaxaDada2$Sum)

TaxaDeblur <- Deblur_table[, c(66, 67)]
colnames(TaxaDeblur)[1] <- "Taxa"
colTaxaDeblur <- ddply(TaxaDeblur, .(Taxa), summarize, Sum=sum(Mean))
colTaxaDeblur$Sum <- colTaxaDeblur$Sum/ sum(colTaxaDeblur$Sum)

TaxaUnoise <- Unoise_table[, c(66, 67)]
colnames(TaxaUnoise)[1] <- "Taxa"

colTaxaUnoise <- ddply(TaxaUnoise, .(Taxa), summarize, Sum=sum(Mean))
colTaxaUnoise$Sum <- colTaxaUnoise$Sum / sum(colTaxaUnoise$Sum)

TaxaOpen <- Open_table[, c(66,67)]
colnames(TaxaOpen)[1] <- "Taxa"
colTaxaOpen <- ddply(TaxaOpen, .(Taxa), summarize, Sum=sum(Mean))
colTaxaOpen$Sum <- colTaxaOpen$Sum / sum(colTaxaOpen$Sum)

colTaxaDada2$group <- "DADA2"
colTaxaDeblur$group <- "Deblur"
colTaxaUnoise$group <- "UNOISE3"
colTaxaOpen$group <- "Open"

colTaxaDada2 <- arrange(colTaxaDada2, desc(Sum))
colTaxaDeblur <- arrange(colTaxaDeblur, desc(Sum))
colTaxaUnoise <- arrange(colTaxaUnoise, desc(Sum))
colTaxaOpen <- arrange(colTaxaOpen, desc(Sum))
colTaxaDada2$order <- seq.int(nrow(colTaxaDada2))
colTaxaDeblur$order <- seq.int(nrow(colTaxaDeblur))
colTaxaUnoise$order <- seq.int(nrow(colTaxaUnoise))
colTaxaOpen$order <- seq.int(nrow(colTaxaOpen))

taxaDat <- rbind(colTaxaDada2, colTaxaDeblur, colTaxaUnoise, colTaxaOpen)

Species <- ggplot(data= taxaDat, aes(x = order, y = log(Sum), color = group)) + geom_line() + labs(y= "log(Relative Abundance)", x = "Rank Species Abundance")+ theme_classic()

grid_plot <- plot_grid(ASVs, Species, labels="AUTO", ncol=2)
grid_plot
