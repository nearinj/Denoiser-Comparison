#make 2d plots for Exercise data

library("ggplot2")
library("vegan")
library("cowplot")
library("gridGraphics")

setwd("/home/jacob/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/fixed_combined/pplacer_distances/")




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
proportions <- read.table("weighted_PC_cords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
SampleNum <- read.table("weighted_PC_cords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
sample_cord <- read.table("weighted_PC_cords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = SampleNum[[2]],
                          header=FALSE, row.names=1)

unifrac_combined_Exercise <- sample_cord[,c("V2", "V3", "V4")]
colnames(unifrac_combined_Exercise) <- c("PC1", "PC2", "PC3")
unifrac_combined_Exercise$sample <- as.factor(gsub("^.+_", "", rownames(unifrac_combined_Exercise)))
unifrac_combined_Exercise$Pipeline <- factor(gsub("_.+$", "", rownames(unifrac_combined_Exercise)))
unifrac_combined_Exercise$soil <- NA
unifrac_combined_Exercise$Pipeline <- gsub("Dada", "DADA2", unifrac_combined_Exercise$Pipeline)
unifrac_combined_Exercise$Pipeline <- gsub("Unoise", "UNOISE3", unifrac_combined_Exercise$Pipeline)
unifrac_combined_Exercise$Pipeline <- gsub("open-ref", "Open", unifrac_combined_Exercise$Pipeline)
unifrac_combined_rm_noPos <- unifrac_combined_Exercise[-which(unifrac_combined_Exercise$Pipeline=="Deblur-noPos"),]

samples_to_keep <- unifrac_combined_rm_noPos[unifrac_combined_rm_noPos$Pipeline=="Deblur", "sample"]

unifrac_combined_filtered <- unifrac_combined_rm_noPos[which(unifrac_combined_rm_noPos$sample %in% samples_to_keep),]
#need to filter the sample list so that there ar enot single samples... 

unifrac_plot <- ggplot(data=unifrac_combined_filtered, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24 ,25)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(unifrac_combined_filtered$sample)) {
  unifrac_plot <- unifrac_plot + geom_line(data = unifrac_combined_filtered[which(unifrac_combined_filtered$sample == sample),], aes(x = PC1, y = PC2))
}

plot(unifrac_plot)


### Plot bray-curtis distances ordination

bray_curtis_dm <- read.table("../../final_combined/bray_curtis/distances/bray_curtis_rare_CombinedV2_Exercise_taxa_L6.txt",
                             sep="\t", 
                             header=T, 
                             stringsAsFactors = FALSE,
                             row.names=1)

#remove Deblur-noPos samples
bray_curtis_dm <- bray_curtis_dm[-(grep("Deblur-noPos.*", rownames(bray_curtis_dm))), ]
bray_curtis_dm <- bray_curtis_dm[,-(grep("Deblur.noPos.*", colnames(bray_curtis_dm)))]

bray_curtis_NMDS <- metaMDS(as.dist(bray_curtis_dm))


bray_curtis_NMDS_df <- data.frame(NMDS1=bray_curtis_NMDS$points[,1], 
                                  NMDS2=bray_curtis_NMDS$points[,2],
                                  sample=as.factor(gsub("^.+_", "", rownames(bray_curtis_NMDS$points))),
                                  Pipeline=as.factor(sub("_.*","",rownames(bray_curtis_NMDS$points))))
bray_curtis_NMDS_df$Pipeline <- gsub("Dada", "DADA2", bray_curtis_NMDS_df$Pipeline)
bray_curtis_NMDS_df$Pipeline <- gsub("open-ref", "Open", bray_curtis_NMDS_df$Pipeline)
bray_curtis_NMDS_df$Pipeline <- gsub("Unoise", "UNOISE3", bray_curtis_NMDS_df$Pipeline)

bray_curtis_NMDS_df <- bray_curtis_NMDS_df[complete.cases(bray_curtis_NMDS_df),]

samples_to_keep <- bray_curtis_NMDS_df[bray_curtis_NMDS_df$Pipeline=="Deblur", "sample"]
bray_curtis_df_filtered <- bray_curtis_NMDS_df[bray_curtis_NMDS_df$sample %in% samples_to_keep,]

#filter out samples removed due to rarefaction


bray_curtis_plot <- ggplot(data=bray_curtis_df_filtered, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)



for (sample in levels(bray_curtis_df_filtered$sample)) {
  bray_curtis_plot <- bray_curtis_plot + geom_line(data = bray_curtis_df_filtered[which(bray_curtis_df_filtered$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

plot(bray_curtis_plot)

#makes figure 4 when other script plots are loaded as well
#par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
#grid_plot <- plot_grid(Unifrac_box, Bray_Curtis_box, unifrac_plot, bray_curtis_plot, labels=c("A", "B", "C", "D"), rel_heights = c(1,2))
#grid_plot

######################################################unweighted unifrac

#read in proportions for the PC
U_proportions <- read.table("unweighted_PC_cords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
U_SampleNum <- read.table("unweighted_PC_cords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
U_sample_cord <- read.table("unweighted_PC_cords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = U_SampleNum[[2]],
                          header=FALSE, row.names=1)

U_unifrac_combined_Exercise <- U_sample_cord[,c("V2", "V3", "V4")]
colnames(U_unifrac_combined_Exercise) <- c("PC1", "PC2", "PC3")
U_unifrac_combined_Exercise$sample <- as.factor(gsub("^.+_", "", rownames(U_unifrac_combined_Exercise)))
U_unifrac_combined_Exercise$Pipeline <- factor(gsub("_.+$", "", rownames(U_unifrac_combined_Exercise)))
U_unifrac_combined_Exercise$Pipeline <- gsub("Dada", "DADA2", U_unifrac_combined_Exercise$Pipeline)
U_unifrac_combined_Exercise$Pipeline <- gsub("Unoise", "UNOISE3", U_unifrac_combined_Exercise$Pipeline)
U_unifrac_combined_Exercise$Pipeline <- gsub("open-ref", "Open", U_unifrac_combined_Exercise$Pipeline)
U_unifrac_combined_rm_noPos <- U_unifrac_combined_Exercise[-which(U_unifrac_combined_Exercise$Pipeline=="Deblur-noPos"),]

U_samples_to_keep <- U_unifrac_combined_rm_noPos[U_unifrac_combined_rm_noPos$Pipeline=="Deblur", "sample"]

U_unifrac_combined_filtered <- U_unifrac_combined_rm_noPos[which(U_unifrac_combined_rm_noPos$sample %in% U_samples_to_keep),]
#need to filter the sample list so that there ar enot single samples... 

U_unifrac_plot <- ggplot(data=U_unifrac_combined_filtered, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(U_proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(U_proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24 ,25)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(U_unifrac_combined_filtered$sample)) {
  U_unifrac_plot <- U_unifrac_plot + geom_line(data = U_unifrac_combined_filtered[which(U_unifrac_combined_filtered$sample == sample),], aes(x = PC1, y = PC2))
}

plot(U_unifrac_plot)
Exercise_U_unifrac_plot <- U_unifrac_plot
