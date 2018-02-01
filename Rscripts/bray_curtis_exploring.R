# Commands used to explore data to try and track down why the Deblur bray-curtis distances were higher.

library(ggplot2)
library(cowplot)
library(gridGraphics)

par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0)) # move plot to the right and up

################## Commands for exploring blueberry results ##################

setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/bray_curtis/")

in_tab <- read.table("tax_sum/merged_blueberry_taxa_L6.txt",
                     header=T,
                     comment.char="",
                     skip=1,
                     sep="\t",
                     stringsAsFactors = FALSE,
                     row.names=1)

# Split by each pipeline.
dada_tab <- in_tab[, grep("Dada", colnames(in_tab))]
deblur_tab <- in_tab[, grep("Deblur", colnames(in_tab))]
unoise_tab <- in_tab[, grep("Unoise", colnames(in_tab))]

dada_tab_means <- rowMeans(dada_tab)
deblur_tab_means <- rowMeans(deblur_tab)
unoise_tab_means <- rowMeans(unoise_tab)

hist(dada_tab_means - deblur_tab_means, xlim=c(-0.03,0.03), breaks=100, main="", xlab="DADA2 - Deblur", ylim=c(0,260), col="grey")
dada_v_deblur_genera_diff <- recordPlot()

hist(dada_tab_means - unoise_tab_means, xlim=c(-0.03,0.03), breaks=50, main="", xlab="DADA2 - UNOISE3", ylim=c(0,260), col="grey")
dada_v_unoise_genera_diff <- recordPlot()

hist(deblur_tab_means - unoise_tab_means, xlim=c(-0.03,0.03), breaks=100, main="", xlab="Deblur - UNOISE3", ylim=c(0,260), col="grey")
deblur_v_unoise_genera_diff <- recordPlot()

plot_grid(dada_v_deblur_genera_diff,
          dada_v_unoise_genera_diff,
          deblur_v_unoise_genera_diff,
          labels="AUTO")

# Identify that taxa that are outliers above abs(0.01) on the DADA2-Deblur and Deblur - UNOISE3 histograms.
dada_deblur_diff <- dada_tab_means - deblur_tab_means
unoise_deblur_diff <- unoise_tab_means - deblur_tab_means

which(abs(dada_deblur_diff) > 0.01)
which(abs(unoise_deblur_diff) > 0.01)

# Plot the 6/7 taxa that are outliers that are above this cut-off for both the DADA2 and UNOISE3 histograms.
boxplot(as.numeric(dada_tab[9,]), 
     as.numeric(deblur_tab[9,]),
     as.numeric(unoise_tab[9,]),
     names = c("Dada", "Deblur", "Unoise"),
     ylab = "Bacteria;Acidobacteria;Acidobacteria_Gp1;Gp1;Unclassified;Unclassified",
     ylim=c(0,0.13))

Gp1_boxplot <- recordPlot()

boxplot(as.numeric(dada_tab[10,]), 
        as.numeric(deblur_tab[10,]),
        as.numeric(unoise_tab[10,]),
        names = c("Dada", "Deblur", "Unoise"),
        ylab = "Bacteria;Acidobacteria;Acidobacteria_Gp1;Granulicella;Unclassified;Unclassified",
        ylim=c(0,0.13))

Granulicella_boxplot <- recordPlot()

# Not going to show this plot, but it also was shared between DADA2/UNOISE3 as an outlier.
# par(mfrow=c(1,1))
# boxplot(as.numeric(dada_tab[271,]), 
#         as.numeric(deblur_tab[271,]),
#         as.numeric(unoise_tab[271,]),
#         names = c("Dada", "Deblur", "Unoise"),
#         ylab = "Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Unclassified;Unclassified ")

boxplot(as.numeric(dada_tab[439,]), 
        as.numeric(deblur_tab[439,]),
        as.numeric(unoise_tab[439,]),
        names = c("Dada", "Deblur", "Unoise"),
        ylab = "Bacteria;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified",
        ylim=c(0,0.23))

Bacteria_Unclassified <- recordPlot()

boxplot(as.numeric(dada_tab[451,]), 
        as.numeric(deblur_tab[451,]),
        as.numeric(unoise_tab[451,]),
        names = c("Dada", "Deblur", "Unoise"),
        ylab = "Bacteria;Verrucomicrobia;Unclassified;Unclassified;Unclassified;Unclassified",
        ylim=c(0,0.23))

Bacteria_Verrucomicrobia <- recordPlot()

boxplot(as.numeric(dada_tab[447,]), 
        as.numeric(deblur_tab[447,]),
        as.numeric(unoise_tab[447,]),
        names = c("Dada", "Deblur", "Unoise"),
        ylab = "Bacteria;Verrucomicrobia;Spartobacteria;Unclassified;Unclassified;Unclassified",
        ylim=c(0,0.23))

Bacteria_Verrucomicrobia_Spartobacteria <- recordPlot()


plot_grid(Bacteria_Unclassified,
          Bacteria_Verrucomicrobia,
          Bacteria_Verrucomicrobia_Spartobacteria,
          Gp1_boxplot,
          Granulicella_boxplot,
          labels="AUTO",
          nrow=2,
          ncol=3)


################## Commands for exploring BISCUIT results ##################

in_tab_biscuit <- read.table("../../../../biscuit/med/COMBINED/Tax_Sum/Combined_biscuit_taxa_L6.txt",
                     header=T,
                     comment.char="",
                     skip=1,
                     sep="\t",
                     stringsAsFactors = FALSE,
                     row.names=1)

# Split by each pipeline.
dada_tab_biscuit <- in_tab_biscuit[, grep("Dada", colnames(in_tab_biscuit))]
deblur_tab_biscuit <- in_tab_biscuit[, grep("Deblur", colnames(in_tab_biscuit))]
unoise_tab_biscuit <- in_tab_biscuit[, grep("Unoise", colnames(in_tab_biscuit))]

dada_tab_biscuit_means <- rowMeans(dada_tab_biscuit)
deblur_tab_biscuit_means <- rowMeans(deblur_tab_biscuit)
unoise_tab_biscuit_means <- rowMeans(unoise_tab_biscuit)

hist(dada_tab_biscuit_means - deblur_tab_biscuit_means, xlim=c(-0.03,0.03), ylim=c(0, 100), breaks=100, main="", xlab="DADA2 - Deblur",  col="grey")
dada_v_deblur_genera_diff_biscuit <- recordPlot()

hist(dada_tab_biscuit_means - unoise_tab_biscuit_means, xlim=c(-0.03,0.03), ylim=c(0, 100), breaks=50, main="", xlab="DADA2 - UNOISE3", col="grey")
dada_v_unoise_genera_diff_biscuit <- recordPlot()

hist(deblur_tab_biscuit_means - unoise_tab_biscuit_means, xlim=c(-0.03,0.03), ylim=c(0, 100), breaks=100, main="", xlab="Deblur - UNOISE3",  col="grey")
deblur_v_unoise_genera_diff_biscuit <- recordPlot()

plot_grid(dada_v_deblur_genera_diff_biscuit,
          dada_v_unoise_genera_diff_biscuit,
          deblur_v_unoise_genera_diff_biscuit,
          labels="AUTO")

dada_deblur_diff_biscuit <- dada_tab_biscuit_means - deblur_tab_biscuit_means
unoise_deblur_diff_biscuit <- unoise_tab_biscuit_means - deblur_tab_biscuit_means
unoise_dada_diff_biscuit <- unoise_tab_biscuit_means - dada_tab_biscuit_means

which(abs(unoise_deblur_diff_biscuit) > 0.005)
which(abs(dada_deblur_diff_biscuit) > 0.005)

### The 2 overlapping both are:
#Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Unclassified 
#72 
#Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia/Shigella 
#146 

boxplot(as.numeric(dada_tab_biscuit[72,]), 
        as.numeric(deblur_tab_biscuit[72,]),
        as.numeric(unoise_tab_biscuit[72,]),
        names = c("Dada", "Deblur", "Unoise"),
        ylab = "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Unclassified",
        ylim=c(0,0.23))

Lachnospiraceae_boxplot <- recordPlot()

boxplot(as.numeric(dada_tab_biscuit[146,]), 
        as.numeric(deblur_tab_biscuit[146,]),
        as.numeric(unoise_tab_biscuit[146,]),
        names = c("Dada", "Deblur", "Unoise"),
        ylab = "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia/Shigella",
        ylim=c(0,0.6))

Ecoli_boxplot <- recordPlot()

plot_grid(Lachnospiraceae_boxplot,
          Ecoli_boxplot,
          labels="AUTO",
          nrow=1,
          ncol=2)

################## The below commented out commands were used for exploring the data, but not for making supplementary figures ################## 

# One sample's rel. abundance in terms of classified and unclassified genera
#all_unclassified <- grep("Unclassified", rownames(in_tab))

# par(mfrow=c(2,3))
# plot(dada_tab[-all_unclassified, "Dada_Bact126"], deblur_tab[-all_unclassified, "Deblur_Bact126"])
# plot(dada_tab[-all_unclassified, "Dada_Bact126"], unoise_tab[-all_unclassified, "Unoise_Bact126"])
# plot(unoise_tab[-all_unclassified, "Unoise_Bact126"], deblur_tab[-all_unclassified, "Deblur_Bact126"])
# 
# plot(dada_tab[all_unclassified, "Dada_Bact126"], deblur_tab[all_unclassified, "Deblur_Bact126"])
# plot(dada_tab[all_unclassified, "Dada_Bact126"], unoise_tab[all_unclassified, "Unoise_Bact126"])
# plot(unoise_tab[all_unclassified, "Unoise_Bact126"], deblur_tab[all_unclassified, "Deblur_Bact126"])
# 
# cor.test(dada_tab[-all_unclassified, "Dada_Bact126"], deblur_tab[-all_unclassified, "Deblur_Bact126"], method="spearman")
# cor.test(dada_tab[all_unclassified, "Dada_Bact126"], deblur_tab[all_unclassified, "Deblur_Bact126"], method="spearman")
# 
# cor.test(dada_tab[-all_unclassified, "Dada_Bact126"], unoise_tab[-all_unclassified, "Unoise_Bact126"], method="spearman")
# cor.test(dada_tab[all_unclassified, "Dada_Bact126"], unoise_tab[all_unclassified, "Unoise_Bact126"], method="spearman")
# 
# cor.test(deblur_tab[-all_unclassified, "Deblur_Bact126"], unoise_tab[-all_unclassified, "Unoise_Bact126"], method="spearman")
# cor.test(deblur_tab[all_unclassified, "Deblur_Bact126"], unoise_tab[all_unclassified, "Unoise_Bact126"], method="spearman")

# # Get row indices for Unclassified taxa at each level and not at above levels.
# k__unclassified <- "Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified"
# in_tab_noK <- in_tab[-which(rownames(in_tab) %in% "Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified"),]
# 
# p__unclassified <- grep("^.+;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified$", rownames(in_tab_noK), value=TRUE)
# in_tab_noK_noP <- in_tab_noK[-which(rownames(in_tab_noK) %in% p__unclassified) , ]
# 
# c__unclassified <- grep("^.+;.+;Unclassified;Unclassified;Unclassified;Unclassified$", rownames(in_tab_noK_noP), value=TRUE)
# in_tab_noK_noP_noC <- in_tab_noK_noP[-which(rownames(in_tab_noK_noP) %in% c__unclassified) , ]
# 
# o__unclassified <- grep("^.+;.+;.+;Unclassified;Unclassified;Unclassified$", rownames(in_tab_noK_noP_noC), value=TRUE)
# in_tab_noK_noP_noC_noO <- in_tab_noK_noP_noC[-which(rownames(in_tab_noK_noP_noC) %in% o__unclassified),]
# 
# f__unclassified <- grep("^.+;.+;.+;.+;Unclassified;Unclassified$", rownames(in_tab_noK_noP_noC_noO), value=TRUE)
# in_tab_noK_noP_noC_noO_noF <- in_tab_noK_noP_noC_noO[-which(rownames(in_tab_noK_noP_noC_noO) %in% f__unclassified),]
# 
# g__unclassified <- grep("^.+;.+;.+;.+;.+;Unclassified$", rownames(in_tab_noK_noP_noC_noO_noF), value=TRUE)
# in_tab_noK_noP_noC_noO_noF_noG <- in_tab_noK_noP_noC_noO_noF[-which(rownames(in_tab_noK_noP_noC_noO_noF) %in% g__unclassified),]
# 

# par(mfrow=c(2,3))
# 
# # Plots of rel. abun. at each unclassified level (means in cases with multiple taxa)
# # Kingdom
# k__unclassified_df <- data.frame(DADA2=as.numeric(dada_tab[k__unclassified,]),
#                                  Deblur=as.numeric(deblur_tab[k__unclassified,]),
#                                  Unoise3=as.numeric(unoise_tab[k__unclassified,]))
# boxplot(k__unclassified_df, ylab = "Rel. Abun. of Unclassified Kingdom", outline=FALSE, ylim=c(0, 0.0028))
# 
# stripchart(k__unclassified_df, add=TRUE,
#            method="jitter", jitter=0.3, vertical = TRUE, pch=1)
# 
# k__dist_plot <- recordPlot()
# 
# kruskal.test(k__unclassified_df)
# #Kruskal-Wallis chi-squared = 4.0431, df = 2, p-value = 0.1324
# 
# # Phylum
# p__unclassified_df <- data.frame(DADA2=as.numeric(dada_tab[p__unclassified,]),
#                                  Deblur=as.numeric(deblur_tab[p__unclassified,]),
#                                  Unoise3=as.numeric(unoise_tab[p__unclassified,]))
# 
# boxplot(p__unclassified_df, ylab = "Rel. Abun. of Unclassified Phylum", outline=FALSE, ylim=c(0, 0.2))
#         
# stripchart(p__unclassified_df, add=TRUE,
#            method="jitter", jitter=0.3, vertical = TRUE, pch=1)
# 
# kruskal.test(p__unclassified_df)
# #Kruskal-Wallis chi-squared = 15.663, df = 2, p-value = 0.000397
# 
# p__dist_plot <- recordPlot()
# 
# # Class
# c__unclassified_df <- data.frame(DADA2=colMeans(dada_tab[c__unclassified,]),
#                                  Deblur=colMeans(deblur_tab[c__unclassified,]),
#                                  Unoise3=colMeans(unoise_tab[c__unclassified,]))
# 
# boxplot(c__unclassified_df, ylab = "Mean Rel. Abun. of Unclassified Classes", outline=FALSE, ylim=c(0, 0.015))
# 
# stripchart(c__unclassified_df, add=TRUE,
#            method="jitter", jitter=0.3, vertical = TRUE, pch=1)
# 
# kruskal.test(c__unclassified_df)
# #Kruskal-Wallis chi-squared = 30.481, df = 2, p-value = 2.405e-07
# 
# c__dist_plot <- recordPlot()
# 
# # Order
# o__unclassified_df <- data.frame(DADA2=colMeans(dada_tab[o__unclassified,]),
#                                  Deblur=colMeans(deblur_tab[o__unclassified,]),
#                                  Unoise3=colMeans(unoise_tab[o__unclassified,]))
# 
# boxplot(o__unclassified_df, ylab = "Mean Rel. Abun. of Unclassified Orders", outline=FALSE, ylim=c(0, 0.02))
# 
# stripchart(o__unclassified_df, add=TRUE,
#            method="jitter", jitter=0.3, vertical = TRUE, pch=1)
# 
# kruskal.test(o__unclassified_df)
# #Kruskal-Wallis chi-squared = 5.7491, df = 2, p-value = 0.05644
# 
# o__dist_plot <- recordPlot()
# 
# # Family
# f__unclassified_df <- data.frame(DADA2=colMeans(dada_tab[f__unclassified,]),
#                                  Deblur=colMeans(deblur_tab[f__unclassified,]),
#                                  Unoise3=colMeans(unoise_tab[f__unclassified,]))
# 
# boxplot(f__unclassified_df, ylab = "Mean Rel. Abun. of Unclassified Families", outline=FALSE, ylim=c(0, 0.0125))
# 
# stripchart(f__unclassified_df, add=TRUE,
#            method="jitter", jitter=0.3, vertical = TRUE, pch=1)
# 
# kruskal.test(f__unclassified_df)
# #Kruskal-Wallis chi-squared = 10.923, df = 2, p-value = 0.004247
# 
# f__dist_plot <- recordPlot()
# 
# # Genus
# g__unclassified_df <- data.frame(DADA2=colMeans(dada_tab[g__unclassified,]),
#                                  Deblur=colMeans(deblur_tab[g__unclassified,]),
#                                  Unoise3=colMeans(unoise_tab[g__unclassified,]))
# 
# boxplot(g__unclassified_df, ylab = "Mean Rel. Abun. of Unclassified Genera", outline=FALSE, ylim=c(0, 0.0025))
# 
# stripchart(g__unclassified_df, add=TRUE,
#            method="jitter", jitter=0.3, vertical = TRUE, pch=1)
# 
# kruskal.test(g__unclassified_df)
# #Kruskal-Wallis chi-squared = 21.254, df = 2, p-value = 2.425e-05
# 
# g__dist_plot <- recordPlot()

# plot_grid(k__dist_plot,
#           p__dist_plot,
#           c__dist_plot,
#           o__dist_plot,
#           f__dist_plot,
#           g__dist_plot,
#           labels=c("A","B","C","D","E","F"))
