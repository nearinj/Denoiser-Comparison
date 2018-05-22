#make 2d plots for blueberry data

library("ggplot2")
library("vegan")
library("cowplot")
library("gridGraphics")

setwd("/media/jacob/Storage/Jacob_Ap_27/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/final_combined/sequences/pynast_aligned_seqs/sep_aligned/")


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


### Plot weighted UniFrac data

#read in proportions for the PC
proportions <- read.table("combined_PC_coords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
SampleNum <- read.table("combined_PC_coords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
sample_cord <- read.table("combined_PC_coords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = SampleNum[[2]],
                          header=FALSE, row.names=1)

unifrac_combined_blueberry <- sample_cord[,c("V2", "V3", "V4")]
colnames(unifrac_combined_blueberry) <- c("PC1", "PC2", "PC3")
unifrac_combined_blueberry$sample <- as.factor(gsub("^.+_", "", rownames(unifrac_combined_blueberry)))
unifrac_combined_blueberry$Pipeline <- factor(gsub("_.+$", "", rownames(unifrac_combined_blueberry)))
unifrac_combined_blueberry$soil <- NA
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% rhizosphere_samples), "soil"] <- "Rhizosphere"
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% bulk_samples), "soil"] <- "Bulk"


#remove noPos samples
#unifrac_combined_blueberry_temp <- unifrac_combined_blueberry[-which(unifrac_combined_blueberry$Pipeline %in% "Deblur-noPos"),]


unifrac_plot <- ggplot(data=unifrac_combined_blueberry, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(unifrac_combined_blueberry$sample)) {
  unifrac_plot <- unifrac_plot + geom_line(data = unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample == sample),], aes(x = PC1, y = PC2))
}

plot(unifrac_plot)


###weighted unifrac original ! with seperatge PC

#read in proportions for the PC
O_proportions <- read.table("weighted_orig_PC_cords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
O_SampleNum <- read.table("weighted_orig_PC_cords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
O_sample_cord <- read.table("weighted_orig_PC_cords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = O_SampleNum[[2]],
                          header=FALSE, row.names=1)

O_unifrac_combined_blueberry <- O_sample_cord[,c("V2", "V3", "V4")]
colnames(O_unifrac_combined_blueberry) <- c("PC1", "PC2", "PC3")
O_unifrac_combined_blueberry$sample <- as.factor(gsub("^.+_", "", rownames(O_unifrac_combined_blueberry)))
O_unifrac_combined_blueberry$Pipeline <- factor(gsub("_.+$", "", rownames(O_unifrac_combined_blueberry)))
O_unifrac_combined_blueberry$soil <- NA
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% rhizosphere_samples), "soil"] <- "Rhizosphere"
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% bulk_samples), "soil"] <- "Bulk"


#remove noPos samples
#unifrac_combined_blueberry_temp <- unifrac_combined_blueberry[-which(unifrac_combined_blueberry$Pipeline %in% "Deblur-noPos"),]


O_unifrac_plot <- ggplot(data=O_unifrac_combined_blueberry, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(O_proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(O_proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(O_unifrac_combined_blueberry$sample)) {
  O_unifrac_plot <- O_unifrac_plot + geom_line(data = O_unifrac_combined_blueberry[which(O_unifrac_combined_blueberry$sample == sample),], aes(x = PC1, y = PC2))
}

plot(O_unifrac_plot)



###### weighted unifrac original

unifrac_combined_blueberry_orig <- unifrac_combined_blueberry[-(which(unifrac_combined_blueberry$Pipeline=="open-ref")),]

unifrac_plot_orig <- ggplot(data=unifrac_combined_blueberry_orig, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)


for (sample in levels(unifrac_combined_blueberry_orig$sample)) {
  unifrac_plot_orig <- unifrac_plot_orig + geom_line(data = unifrac_combined_blueberry_orig[which(unifrac_combined_blueberry_orig$sample == sample),], aes(x = PC1, y = PC2))
}

unifrac_plot_orig

######################################################unweighted plot
#read in proportions for the PC
U_proportions <- read.table("unweighted_PC_cords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
U_SampleNum <- read.table("unweighted_PC_cords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
U_sample_cord <- read.table("unweighted_PC_cords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = U_SampleNum[[2]],
                          header=FALSE, row.names=1)

U_unifrac_combined_blueberry <- U_sample_cord[,c("V2", "V3", "V4")]
colnames(U_unifrac_combined_blueberry) <- c("PC1", "PC2", "PC3")
U_unifrac_combined_blueberry$sample <- as.factor(gsub("^.+_", "", rownames(U_unifrac_combined_blueberry)))
U_unifrac_combined_blueberry$Pipeline <- factor(gsub("_.+$", "", rownames(U_unifrac_combined_blueberry)))
U_unifrac_combined_blueberry$soil <- NA
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% rhizosphere_samples), "soil"] <- "Rhizosphere"
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% bulk_samples), "soil"] <- "Bulk"


#remove noPos samples
#unifrac_combined_blueberry_temp <- unifrac_combined_blueberry[-which(unifrac_combined_blueberry$Pipeline %in% "Deblur-noPos"),]


U_unifrac_plot <- ggplot(data=U_unifrac_combined_blueberry, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(U_proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(U_proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(U_unifrac_combined_blueberry$sample)) {
  U_unifrac_plot <- U_unifrac_plot + geom_line(data = U_unifrac_combined_blueberry[which(U_unifrac_combined_blueberry$sample == sample),], aes(x = PC1, y = PC2))
}

plot(U_unifrac_plot)

###################### unweighted original

U_unifrac_combined_blueberry_orig <- U_unifrac_combined_blueberry[-(which(U_unifrac_combined_blueberry$Pipeline=="open-ref")),]

U_unifrac_plot_orig <- ggplot(data=U_unifrac_combined_blueberry_orig, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(U_proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(U_proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)


for (sample in levels(U_unifrac_combined_blueberry_orig$sample)) {
  U_unifrac_plot_orig <- U_unifrac_plot_orig + geom_line(data = U_unifrac_combined_blueberry_orig[which(U_unifrac_combined_blueberry_orig$sample == sample),], aes(x = PC1, y = PC2))
}

U_unifrac_plot_orig

### Plot bray-curtis distances ordination

bray_curtis_dm <- read.table("../../../taxa_summary/bray_curtis/bray_curtis_rare_CombinedV2_Blueberry_taxa_L6.txt",
                             sep="\t", 
                             header=T, 
                             stringsAsFactors = FALSE,
                             row.names=1)
colnames(bray_curtis_dm) <- gsub("001101", "", colnames(bray_curtis_dm))
rownames(bray_curtis_dm) <- gsub("001101", "", rownames(bray_curtis_dm))


#bray_curtis_NMDS <- metaMDS(as.dist(bray_curtis_dm), trymax=10000, k=2, maxit=100)
#saveRDS(bray_curtis_NMDS, file="bray_curtis_NMDs.rds")
bray_curtis_NMDS <- readRDS("bray_curtis_NMDs.rds")



bray_curtis_NMDS_df <- data.frame(NMDS1=bray_curtis_NMDS$points[,1], 
                                 NMDS2=bray_curtis_NMDS$points[,2],
                                 sample=as.factor(gsub("^.+_", "", rownames(bray_curtis_NMDS$points))),
                                 Pipeline=as.factor(sub("_.*", "", rownames(bray_curtis_NMDS$points))))

bray_curtis_NMDS_df <- bray_curtis_NMDS_df[complete.cases(bray_curtis_NMDS_df),]

bray_curtis_NMDS_df_NO_noPos <- bray_curtis_NMDS_df[-which((bray_curtis_NMDS_df$Pipeline == "Deblur-noPos")),]

bray_curtis_plot <- ggplot(data=bray_curtis_NMDS_df_NO_noPos, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)



for (sample in levels(bray_curtis_NMDS_df_NO_noPos$sample)) {
  bray_curtis_plot <- bray_curtis_plot + geom_line(data = bray_curtis_NMDS_df_NO_noPos[which(bray_curtis_NMDS_df_NO_noPos$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

plot(bray_curtis_plot)

### bray curtis original

bray_curtis_NMDS_df_orig <- bray_curtis_NMDS_df_NO_noPos[-(which(bray_curtis_NMDS_df_NO_noPos$Pipeline=="open-ref")),]

bray_curtis_plot_orig <- ggplot(data=bray_curtis_NMDS_df_orig, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(bray_curtis_NMDS_df_orig$sample)) {
  bray_curtis_plot_orig <- bray_curtis_plot_orig + geom_line(data = bray_curtis_NMDS_df_orig[which(bray_curtis_NMDS_df_orig$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

bray_curtis_plot_orig


#makes figure 4 when other script plots are loaded as well
#par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
#grid_plot <- plot_grid(Unifrac_box, Bray_Curtis_box, unifrac_plot, bray_curtis_plot, labels=c("A", "B", "C", "D"))
#grid_plot

