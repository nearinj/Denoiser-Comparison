#make 2d plots for blueberry

setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/")
#read in proportions for the PC
proportions <- read.table("plots/bdiv-out/weighted_unifrac_pc.txt", sep = "\t", nrow=1, skip=4)
#read in the number oof samples
SampleNum <- read.table("plots/bdiv-out/weighted_unifrac_pc.txt", sep = "\t", nrow=1)

#set cords for samplles
blueberry_meta <- read.table("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/map_blueberry_merged.csv",
                             header=T,
                             comment.char="",
                             stringsAsFactors = FALSE,
                             sep="\t")

rhizosphere_samples <- unique(gsub(".*_", "", blueberry_meta[which(blueberry_meta$Description=="Rhizosphere"), "SampleID"]))
bulk_samples <- unique(gsub(".*_", "", blueberry_meta[which(blueberry_meta$Description=="Bulk"), "SampleID"]))

sample_cord <- read.table("plots/bdiv-out/weighted_unifrac_pc.txt", 
                               sep ="\t", skip = 9, 
                               nrow = SampleNum[[2]],
                               header=FALSE, row.names=1)

Shapes <- rep(x = 0, nrow(sample_cord))
Colors <- rep(x = "green", nrow(sample_cord))

Shapes[grep("Dada", rownames(sample_cord)) ] <- 15
Shapes[grep("Unoise", rownames(sample_cord))] <- 16
Shapes[grep("Deblur", rownames(sample_cord))] <- 17

for (val in Rhizosphere){
  Colors[grep(val, rownames(sample_cord))] <- "red"
}

paste(Colors)
plot(sample_cord$V2, sample_cord$V3, 
     xlab=paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep=""),  
     ylab=paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep=""),
     pch=c(Shapes), bg=Colors, col=Colors)

# Colours taken from here: http://godsnotwheregodsnot.blogspot.ca/2012/09/color-distribution-methodology.html
diff_col <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
              "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
              "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
              "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
              "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
              "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
              "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
              "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
              "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
              "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
              "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
              "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

library("ggplot2")

unifrac_combined_blueberry <- sample_cord[,c("V2", "V3", "V4")]
colnames(unifrac_combined_blueberry) <- c("PC1", "PC2", "PC3")
unifrac_combined_blueberry$sample <- as.factor(gsub("^.+_", "", rownames(unifrac_combined_blueberry)))
unifrac_combined_blueberry$Pipeline <- factor(gsub("_.+$", "", rownames(unifrac_combined_blueberry)))
unifrac_combined_blueberry$soil <- NA
unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% rhizosphere_samples), "soil"] <- "Rhizosphere"
unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% bulk_samples), "soil"] <- "Bulk"

#unifrac_combined_blueberry$sample_numeric <- as.numeric(unifrac_combined_blueberry$sample)*100

p <- ggplot(data=unifrac_combined_blueberry, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(unifrac_combined_blueberry$sample)) {
  p <- p + geom_line(data = unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample == sample),], aes(x = PC1, y = PC2))
}

plot(p)



library("vegan")
bray_curtis_dm <- read.table("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/bray_curtis/tax_sum/beta_diversity/bray_curtis_merged_blueberry_taxa_L6.txt",
                             sep="\t", 
                             header=T, 
                             stringsAsFactors = FALSE,
                             row.names=1)

bray_curtis_NMDS <- metaMDS(as.dist(bray_curtis_dm))
rownames(blueberry_meta) <- blueberry_meta$SampleID

bray_curtis_NMDS_df <- data.frame(NMDS1=bray_curtis_NMDS$points[,1], 
                                 NMDS2=bray_curtis_NMDS$points[,2],
                                 sample=as.factor(gsub("^.+_", "", rownames(bray_curtis_NMDS$points))),
                                 Pipeline=as.factor(blueberry_meta[rownames(bray_curtis_NMDS$points), "Pipeline"]))

bray_curtis_plot <- ggplot(data=bray_curtis_NMDS_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(bray_curtis_NMDS_df$sample)) {
  bray_curtis_plot <- bray_curtis_plot + geom_line(data = bray_curtis_NMDS_df[which(bray_curtis_NMDS_df$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

plot(bray_curtis_plot)



combined_genus <- read.table("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/bray_curtis/tax_sum/merged_blueberry_taxa_L6.txt",
                             sep="\t", 
                             header=T, 
                             stringsAsFactors = FALSE,
                             skip=1,
                             comment.char="",
                             row.names=1)

combined_genus_nonzero <- colSums(combined_genus > 0)

boxplot(combined_genus_nonzero[grep("Deblur", names(combined_genus_nonzero))],
        combined_genus_nonzero[grep("Unoise", names(combined_genus_nonzero))],
        combined_genus_nonzero[grep("Dada", names(combined_genus_nonzero))],
        names=c("Deblur", "Unoise", "Dada"))