#Bray-Curtis figure with all three datasets on it!

library(ggplot2)

### Biscuit boxplot!
setwd("/home/jacob/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/final_combined/bray_curtis/distances/")

par=(mfrow=c(2,2))
# Read in bray-curtis distance matrix (at genus level) of samples from all 3 pipelines
Bi_DM <- read.table("bray_curtis_rare_CombinedV2_biscuit_taxa_L6.txt", sep="\t", header=TRUE, row.names=1, check.names = F)

# Read in mapfile.
#map <- read.table("../../map.txt", sep="\t", header=T, stringsAsFactors = FALSE,) 

# Function that will reorder rows and columns to be the same order.
reorder_pipeline_samples <- function(dm_input, row_prefix, col_prefix) {
  
  dm_input_rownames_sample <- gsub(row_prefix, "", rownames(dm_input))
  colnames(dm_input) <- gsub(col_prefix, "", colnames(dm_input))
  dm_input_ordered <- dm_input[, dm_input_rownames_sample]
  colnames(dm_input_ordered) <- gsub("^", col_prefix, colnames(dm_input_ordered))
  
  dm_input_colnames_sample <- gsub(col_prefix, "", colnames(dm_input_ordered))
  
  if(! identical(dm_input_colnames_sample, dm_input_rownames_sample)) {
    stop("Error final row and column names don't match")
  } else {
    return(dm_input_ordered)
  }
}

# Unoise vs Dada
Bi_unoiseVDada <- Bi_DM[grepl("Unoise*", names(Bi_DM)) ,]
Bi_unoiseVDada <- Bi_unoiseVDada[, grepl("Dada*", names(Bi_unoiseVDada))]
Bi_UvDa <- data.matrix(Bi_unoiseVDada)
# Already in order so don't need to re-order this matrix.

# Unoise vs Deblur
Bi_unoiseVDeblur <- Bi_DM[grepl("Unoise*", names(Bi_DM)),]
Bi_unoiseVDeblur <- Bi_unoiseVDeblur[, grepl("Deblur*", names(Bi_unoiseVDeblur))]
Bi_UvDe <- data.matrix(Bi_unoiseVDeblur)

Bi_UvDe_ordered <- reorder_pipeline_samples(dm_input=Bi_UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

# Unoise vs Deblur
Bi_dadaVdeblur <- Bi_DM[grepl("Dada*", names(Bi_DM)) ,]
Bi_dadaVdeblur <- Bi_dadaVdeblur[,grepl("Deblur*", names(Bi_dadaVdeblur))]
Bi_DavDe <- data.matrix(Bi_dadaVdeblur)

Bi_DavDe_ordered <- reorder_pipeline_samples(dm_input=Bi_DavDe, row_prefix="Dada_", col_prefix="Deblur_")

#DeblurVnoPos
Bi_DeblurVnoPos <- Bi_DM[grep("Deblur_", names(Bi_DM)), grep("Deblur-noPos_", names(Bi_DM))]
colnames(Bi_DeblurVnoPos) <- gsub("001101", "", colnames(Bi_DeblurVnoPos))
Bi_DeVnoPos <- data.matrix(Bi_DeblurVnoPos)
Bi_DeVnoPos_ordered <- reorder_pipeline_samples(dm_input = Bi_DeVnoPos, row_prefix = "Deblur_", col_prefix = "Deblur-noPos_")


#DeblurVOpen
Bi_DeblurVOpen <- data.matrix(Bi_DM[grep("Deblur_", names(Bi_DM)), grep("open-ref_", names(Bi_DM))])
Bi_DeVOp <- reorder_pipeline_samples(Bi_DeblurVOpen, row_prefix = "Deblur_", col_prefix = "open-ref_")


#DadaVOpen
Bi_DadaVOpen <- data.matrix(Bi_DM[grep("Dada_", names(Bi_DM)), grep("open-ref_", names(Bi_DM))])
Bi_DaVOp <- reorder_pipeline_samples(Bi_DadaVOpen, row_prefix = "Dada_", col_prefix = "open-ref_")

#UnoiseVOpen
Bi_UnoiseVOpen <- data.matrix(Bi_DM[grep("Unoise_", names(Bi_DM)), grep("open-ref_", names(Bi_DM))])
Bi_UvOP <- reorder_pipeline_samples(Bi_UnoiseVOpen, row_prefix = "Unoise__", col_prefix = "open-ref_")

# Sanity checks:        diag(UnoiseVOpen_order),
dim(Bi_UvDa)
dim(Bi_UvDe_ordered)
dim(Bi_DavDe_ordered)
dim(Bi_DeVnoPos_ordered)
dim(Bi_DeVOp)
dim(Bi_DaVOp)
dim(Bi_UvOP)

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Bi_Bray_Curits_Boxplot <- boxplot(diag(Bi_DavDe_ordered),
                               diag(Bi_UvDa),
                               diag(Bi_UvDe_ordered),
                               diag(Bi_DaVOp),
                               diag(Bi_DeVOp),
                               diag(Bi_UvOP),
                               names=c("DADA2<->Deblur", "DADA2<->UNOISE3", "UNOISE3<->Deblur", "DADA2<->OTU", "Deblur<->OTU", "UNOISE3<->OTU"),
                               ylab="Bray-Curtis Distance",
                               cex.axis=0.66,
                               outline=FALSE)

#Used to make figure 4
Bi_Bray_Curtis_box <- recordPlot(Bi_Bray_Curits_Boxplot)

#### ordination Biscuit
library(vegan)

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

Bi_bray_curtis_NMDS <- metaMDS(as.dist(Bi_DM))


Bi_bray_curtis_NMDS_df <- data.frame(NMDS1=Bi_bray_curtis_NMDS$points[,1], 
                                  NMDS2=Bi_bray_curtis_NMDS$points[,2],
                                  sample=as.factor(gsub("^.+_", "", rownames(Bi_bray_curtis_NMDS$points))),
                                  Pipeline=as.factor(sub("_.*","",rownames(Bi_bray_curtis_NMDS$points))))

Bi_bray_curtis_NMDS_df <- Bi_bray_curtis_NMDS_df[complete.cases(Bi_bray_curtis_NMDS_df),]
Bi_bray_curtis_NMDS_df$Pipeline <- gsub("Unoise", "UNOISE3", Bi_bray_curtis_NMDS_df$Pipeline)
Bi_bray_curtis_NMDS_df$Pipeline <- gsub("Dada", "DADA2", Bi_bray_curtis_NMDS_df$Pipeline)
Bi_bray_curtis_NMDS_df$Pipeline <- gsub("open-ref", "OTU", Bi_bray_curtis_NMDS_df$Pipeline)

#remove noPos
Bi_bray_curtis_NMDS_df_NO_noPos <- Bi_bray_curtis_NMDS_df[-which(Bi_bray_curtis_NMDS_df$Pipeline == "Deblur-noPos"),]
#remove samples that are not in all the others after rarefaction
Bi_samples_to_keep <- Bi_bray_curtis_NMDS_df_NO_noPos[which(Bi_bray_curtis_NMDS_df$Pipeline == "Deblur"), "sample"]

Bi_bray_curtis_NMDS_df_NO_noPos <- Bi_bray_curtis_NMDS_df_NO_noPos[Bi_bray_curtis_NMDS_df_NO_noPos$sample %in% as.character(Bi_samples_to_keep),]

Bi_bray_curtis_plot <- ggplot(data=Bi_bray_curtis_NMDS_df_NO_noPos, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)



for (sample in levels(Bi_bray_curtis_NMDS_df_NO_noPos$sample)) {
  Bi_bray_curtis_plot <- Bi_bray_curtis_plot + geom_line(data = Bi_bray_curtis_NMDS_df_NO_noPos[which(Bi_bray_curtis_NMDS_df_NO_noPos$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

plot(Bi_bray_curtis_plot)




##### Exercise boxplot Bray

setwd("/home/jacob/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/final_combined/bray_curtis/distances/")

par=(mfrow=c(2,2))
# Read in bray-curtis distance matrix (at genus level) of samples from all 3 pipelines
E_DM <- read.table("bray_curtis_rare_CombinedV2_Exercise_taxa_L6.txt", sep="\t", header=TRUE, row.names=1, check.names=F)


# Unoise vs Dada
E_unoiseVDada <- E_DM[grepl("Unoise*", names(E_DM)) ,]

E_unoiseVDada <- E_unoiseVDada[, grepl("Dada*", names(E_unoiseVDada))]
E_UvDa <- data.matrix(E_unoiseVDada)
# Already in order so don't need to re-order this matrix.

# Unoise vs Deblur
E_unoiseVDeblur <- E_DM[grepl("Unoise*", names(E_DM)),]
E_unoiseVDeblur <- E_unoiseVDeblur[, grepl("Deblur*", names(E_unoiseVDeblur))]
E_UvDe <- data.matrix(E_unoiseVDeblur)

E_UvDe_ordered <- reorder_pipeline_samples(dm_input=E_UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

# Unoise vs Deblur
E_dadaVdeblur <- E_DM[grepl("Dada*", names(E_DM)) ,]
E_dadaVdeblur <- E_dadaVdeblur[,grepl("Deblur*", names(E_dadaVdeblur))]
E_DavDe <- data.matrix(E_dadaVdeblur)

E_DavDe_ordered <- reorder_pipeline_samples(dm_input=E_DavDe, row_prefix="Dada_", col_prefix="Deblur_")


# Sanity checks:
dim(E_UvDa)
dim(E_UvDe_ordered)
dim(E_DavDe_ordered)

#DeblurVnoPos
E_DeblurVnoPos <- as.matrix(E_DM[grep("Deblur_", names(E_DM)), grep("Deblur-noPos_", names(E_DM))])
colnames(E_DeblurVnoPos) <- gsub("001101", "", colnames(E_DeblurVnoPos))
E_DeVnoPos <- reorder_pipeline_samples(E_DeblurVnoPos, "Deblur_", "Deblur-noPos_")

#Deblur V Open
E_DeblurVOpen <- as.matrix(E_DM[grep("Deblur_", names(E_DM)), grep("open-ref_", names(E_DM))])
E_DeVOP <- reorder_pipeline_samples(E_DeblurVOpen, "Deblur_", "open-ref_")

#Dada V Open
E_DadaVOpen <- as.matrix(E_DM[grep("Dada_", names(E_DM)), grep("open-ref_", names(E_DM))])
E_DaVOP <- reorder_pipeline_samples(E_DadaVOpen, "Dada_", "open-ref_")

#Unoise V Open

E_UnoiseVOpen <- as.matrix(E_DM[grep("Unoise_", names(E_DM)), grep("open-ref_", names(E_DM))])
E_UVOP <- reorder_pipeline_samples(E_UnoiseVOpen, "Unoise__", "open-ref_")

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
E_Bray_Curits_Boxplot <- boxplot(diag(E_DavDe_ordered),
                               diag(E_UvDa),
                               diag(E_UvDe_ordered),
                               diag(E_DaVOP),
                               diag(E_DeVOP),
                               diag(E_UVOP),
                               names=c("DADA2<->Deblur", "DADA2<->UNOISE3", "UNOISE3<->Deblur", "DADA2<->OTU", "Debur<->OTU", "UNOISE3<->OTU"),
                               ylab="Bray-Curtis Distance",
                               cex.axis=0.66,
                               outline=FALSE)

#Used to make figure 4
E_Bray_Curtis_box <- recordPlot(E_Bray_Curits_Boxplot)


##### Exercise ordination

E_bray_curtis_dm <- read.table("bray_curtis_rare_CombinedV2_Exercise_taxa_L6.txt",
                             sep="\t", 
                             header=T, 
                             stringsAsFactors = FALSE,
                             row.names=1)

#remove Deblur-noPos samples
E_bray_curtis_dm <- E_bray_curtis_dm[-(grep("Deblur-noPos.*", rownames(E_bray_curtis_dm))), ]
E_bray_curtis_dm <- E_bray_curtis_dm[,-(grep("Deblur.noPos.*", colnames(E_bray_curtis_dm)))]

E_bray_curtis_NMDS <- metaMDS(as.dist(E_bray_curtis_dm))


E_bray_curtis_NMDS_df <- data.frame(NMDS1=E_bray_curtis_NMDS$points[,1], 
                                  NMDS2=E_bray_curtis_NMDS$points[,2],
                                  sample=as.factor(gsub("^.+_", "", rownames(E_bray_curtis_NMDS$points))),
                                  Pipeline=as.factor(sub("_.*","",rownames(E_bray_curtis_NMDS$points))))

E_bray_curtis_NMDS_df$Pipeline <- gsub("Dada", "DADA2", E_bray_curtis_NMDS_df$Pipeline)
E_bray_curtis_NMDS_df$Pipeline <- gsub("open-ref", "OTU", E_bray_curtis_NMDS_df$Pipeline)
E_bray_curtis_NMDS_df$Pipeline <- gsub("Unoise", "UNOISE3", E_bray_curtis_NMDS_df$Pipeline)

E_bray_curtis_NMDS_df <- E_bray_curtis_NMDS_df[complete.cases(E_bray_curtis_NMDS_df),]

E_samples_to_keep <- E_bray_curtis_NMDS_df[E_bray_curtis_NMDS_df$Pipeline=="Deblur", "sample"]
E_bray_curtis_df_filtered <- E_bray_curtis_NMDS_df[E_bray_curtis_NMDS_df$sample %in% E_samples_to_keep,]

#filter out samples removed due to rarefaction


E_bray_curtis_plot <- ggplot(data=E_bray_curtis_df_filtered, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)



for (sample in levels(E_bray_curtis_df_filtered$sample)) {
  E_bray_curtis_plot <- E_bray_curtis_plot + geom_line(data = E_bray_curtis_df_filtered[which(E_bray_curtis_df_filtered$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

plot(E_bray_curtis_plot)




#### Blueberry bray-curtis boxplot

setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/final_combined/taxa_summary/bray_curtis/")

par=(mfrow=c(2,2))
# Read in bray-curtis distance matrix (at genus level) of samples from all 3 pipelines
BL_DM <- read.table("bray_curtis_rare_CombinedV2_Blueberry_taxa_L6.txt", sep="\t", header=TRUE, row.names=1, check.names = F)

# Read in mapfile.
#bluberry_map <- read.table("../map_blueberry_merged.csv", sep="\t", header=T, stringsAsFactors = FALSE,) 

# Function that will reorder rows and columns to be the same order.
reorder_pipeline_samples <- function(dm_input, row_prefix, col_prefix) {
  
  dm_input_rownames_sample <- gsub(row_prefix, "", rownames(dm_input))
  colnames(dm_input) <- gsub(col_prefix, "", colnames(dm_input))
  dm_input_ordered <- dm_input[, dm_input_rownames_sample]
  colnames(dm_input_ordered) <- gsub("^", col_prefix, colnames(dm_input_ordered))
  
  dm_input_colnames_sample <- gsub(col_prefix, "", colnames(dm_input_ordered))
  
  if(! identical(dm_input_colnames_sample, dm_input_rownames_sample)) {
    stop("Error final row and column names don't match")
  } else {
    return(dm_input_ordered)
  }
}

# Unoise vs Dada
BL_unoiseVDada <- BL_DM[grepl("Unoise*", names(BL_DM)) ,]
BL_unoiseVDada <- BL_unoiseVDada[, grepl("Dada*", names(BL_unoiseVDada))]
BL_UvDa <- data.matrix(BL_unoiseVDada)
# Already in order so don't need to re-order this matrix.

# Unoise vs Deblur
BL_unoiseVDeblur <- BL_DM[grepl("Unoise*", names(BL_DM)),]
BL_unoiseVDeblur <- BL_unoiseVDeblur[, grepl("Deblur*", names(BL_unoiseVDeblur))]
BL_UvDe <- data.matrix(BL_unoiseVDeblur)

BL_UvDe_ordered <- reorder_pipeline_samples(dm_input=BL_UvDe, row_prefix="Unoise_", col_prefix="Deblur_")

# Unoise vs Deblur
BL_dadaVdeblur <- BL_DM[grepl("Dada*", names(BL_DM)) ,]
BL_dadaVdeblur <- BL_dadaVdeblur[,grepl("Deblur*", names(BL_dadaVdeblur))]
BL_DavDe <- data.matrix(BL_dadaVdeblur)

BL_DavDe_ordered <- reorder_pipeline_samples(dm_input=BL_DavDe, row_prefix="Dada_", col_prefix="Deblur_")


# Sanity checks:
dim(BL_UvDa)
dim(BL_UvDe_ordered)
dim(BL_DavDe_ordered)

BL_DM["Dada_Bact194", "Deblur_Bact194"]
BL_DavDe_ordered["Dada_Bact194", "Deblur_Bact194"]
BL_DM["Unoise_Bact1z2z2","Dada_Bact1z2z2"]
BL_UvDa["Unoise_Bact1z2z2",which(rownames(BL_UvDa) == "Unoise_Bact1z2z2")]

BL_DM["Unoise_Bact1z2z2","Deblur_Bact1z2z2"]
BL_UvDe_ordered["Unoise_Bact1z2z2",which(rownames(BL_UvDe_ordered) == "Unoise_Bact1z2z2")]

#DeblurVnoPos
BL_DeblurVnoPos <- as.matrix(BL_DM[grep("Deblur_", names(BL_DM)), grep("Deblur-noPos_", names(BL_DM))])
colnames(BL_DeblurVnoPos) <- gsub("001101", "", colnames(BL_DeblurVnoPos))
BL_DeVnoPos <- reorder_pipeline_samples(BL_DeblurVnoPos, row_prefix = "Deblur_", col_prefix = "Deblur-noPos_")

#DeblurVOpen
BL_DeblurVOpen <- as.matrix(BL_DM[grep("Deblur_", names(BL_DM)), grep("open-ref_", names(BL_DM))])
BL_DeVOP <- reorder_pipeline_samples(BL_DeblurVOpen, "Deblur_", "open-ref_")

#DadaVOpen
BL_DadaVOpen <- as.matrix(BL_DM[grep("Dada_", names(BL_DM)), grep("open-ref_", names(BL_DM))])
BL_DaVOP <- reorder_pipeline_samples(BL_DadaVOpen, "Dada_", "open-ref_")

#UnoiseVOpen
BL_UnoiseVOpen <- as.matrix(BL_DM[grep("Unoise_", names(BL_DM)), grep("open-ref_", names(BL_DM))])
BL_UVOP <- reorder_pipeline_samples(BL_UnoiseVOpen, "Unoise_", "open-ref_")

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
BL_Bray_Curits_Boxplot <- boxplot(diag(BL_DavDe_ordered),
                               diag(BL_UvDa),
                               diag(BL_UvDe_ordered),
                               diag(BL_DaVOP),
                               diag(BL_DeVOP),
                               diag(BL_UVOP),
                               names=c("DADA2<->Deblur", "DADA2<->UNOISE3", "UNOISE3<->Deblur", "DADA2<->OTU", "Deblur<->OTU", "UNOISE3<->OTU"),
                               ylab="Bray-Curtis Distance", 
                               cex.axis=0.66,
                               outline=FALSE)

#Used to make figure 4
BL_Bray_Curtis_box <- recordPlot(BL_Bray_Curits_Boxplot)




### Blueberry Bray-curtis ordination


BL_bray_curtis_NMDS <- readRDS("../../sequences/pynast_aligned_seqs/sep_aligned/bray_curtis_NMDs.rds")



BL_bray_curtis_NMDS_df <- data.frame(NMDS1=BL_bray_curtis_NMDS$points[,1], 
                                  NMDS2=BL_bray_curtis_NMDS$points[,2],
                                  sample=as.factor(gsub("^.+_", "", rownames(BL_bray_curtis_NMDS$points))),
                                  Pipeline=as.factor(sub("_.*", "", rownames(BL_bray_curtis_NMDS$points))))

BL_bray_curtis_NMDS_df <- BL_bray_curtis_NMDS_df[complete.cases(BL_bray_curtis_NMDS_df),]

BL_bray_curtis_NMDS_df_NO_noPos <- BL_bray_curtis_NMDS_df[-which((BL_bray_curtis_NMDS_df$Pipeline == "Deblur-noPos")),]

BL_bray_curtis_NMDS_df_NO_noPos$Pipeline <- gsub("Dada", "DADA2", BL_bray_curtis_NMDS_df_NO_noPos$Pipeline)
BL_bray_curtis_NMDS_df_NO_noPos$Pipeline <- gsub("open-ref", "OTU", BL_bray_curtis_NMDS_df_NO_noPos$Pipeline)
BL_bray_curtis_NMDS_df_NO_noPos$Pipeline <- gsub("Unoise", "UNOISE3", BL_bray_curtis_NMDS_df_NO_noPos$Pipeline)

BL_bray_curtis_plot <- ggplot(data=BL_bray_curtis_NMDS_df_NO_noPos, aes(NMDS1, NMDS2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)



for (sample in levels(BL_bray_curtis_NMDS_df_NO_noPos$sample)) {
  BL_bray_curtis_plot <- BL_bray_curtis_plot + geom_line(data = BL_bray_curtis_NMDS_df_NO_noPos[which(BL_bray_curtis_NMDS_df_NO_noPos$sample == sample),], aes(x = NMDS1, y = NMDS2))
}

plot(BL_bray_curtis_plot)




library(cowplot)

bray_curtis_fig <- plot_grid(BL_Bray_Curtis_box, Bi_Bray_Curtis_box, E_Bray_Curtis_box, BL_bray_curtis_plot, Bi_bray_curtis_plot, E_bray_curtis_plot,
                             labels="AUTO", rel_heights = c(1,2))
bray_curtis_fig


