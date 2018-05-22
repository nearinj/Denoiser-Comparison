#New Figure 3 with all the weighted Unifrac Values on it!

#Load in blueberry weighted info.

library(ggplot2)

setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/fixed_combined/pplacer_distances/")
par=(mfrow=c(2,2))

# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
DM <- read.table("weighted_unifrac_rare_Combined_Blueberry_fixed.txt", sep="\t", header=TRUE, row.names=1)
rownames(DM) <- gsub("-", ".", rownames(DM))


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
unoiseVDada <- DM[grepl("Unoise*", names(DM)) ,]
unoiseVDada <- unoiseVDada[, grepl("Dada*", names(unoiseVDada))]
UvDa <- data.matrix(unoiseVDada)
# Already in order so don't need to re-order!

# Unoise vs Deblur
unoiseVDeblur <- DM[grepl("Unoise*", names(DM)),]
unoiseVDeblur <- unoiseVDeblur[, grepl("Deblur*", names(unoiseVDeblur))]
UvDe <- data.matrix(unoiseVDeblur)

UvDe_ordered <- reorder_pipeline_samples(dm_input=UvDe, row_prefix="Unoise_", col_prefix="Deblur_")

# Unoise vs Deblur
dadaVdeblur <- DM[grepl("Dada*", names(DM)) ,]
dadaVdeblur <- dadaVdeblur[,grepl("Deblur*", names(dadaVdeblur))]
DavDe <- data.matrix(dadaVdeblur)

DavDe_ordered <- reorder_pipeline_samples(dm_input=DavDe, row_prefix="Dada_", col_prefix="Deblur_")

# Dada vs Dada
dadaVdada <- DM[grep("Dada", names(DM)) ,grep("Dada", names(DM))]
DavDa <- data.matrix(dadaVdada)
DavDa_ordered <- reorder_pipeline_samples(dm_input=DavDa, row_prefix="Dada_", col_prefix="Dada_")

# Deblur vs Deblur
deblurVdeblur <- DM[grep("Deblur_", names(DM)) ,grep("Deblur_", names(DM))]
DevDe <- data.matrix(deblurVdeblur)
DevDe_ordered <- reorder_pipeline_samples(dm_input=DevDe, row_prefix="Deblur_", col_prefix="Deblur_")

# unoise vs unoise
unoiseVunoise <- DM[grep("Unoise", names(DM)) ,grep("Unoise", names(DM))]
UvU <- data.matrix(unoiseVunoise)
UvU_ordered <- reorder_pipeline_samples(dm_input=UvU, row_prefix="Unoise_", col_prefix="Unoise_")


#add in open
DadaVOpen <- DM[grep("open.ref_", names(DM)), grep("Dada_", names(DM))]
DadaVOpenorder <- reorder_pipeline_samples(dm_input=as.matrix(DadaVOpen), row_prefix="open.ref_", col_prefix ="Dada_" )

DeblurVOpen <- DM[grep("open.ref_", names(DM)), grep("Deblur_", names(DM))]
DeblurVOpen_order <- reorder_pipeline_samples(dm_input = as.matrix(DeblurVOpen), row_prefix = "open.ref_", col_prefix = "Deblur_")

UnoiseVOpen <- DM[grep("open.ref_", names(DM)), grep("Unoise_", names(DM))]
UnoiseVOpen_order <- reorder_pipeline_samples(dm_input = as.matrix(UnoiseVOpen), row_prefix = "open.ref_", col_prefix = "Unoise_")



# Sanity checks:
dim(UvDa)
dim(UvDe_ordered)
dim(DavDe_ordered)
dim(DadaVOpenorder)
dim(DeblurVOpen_order)
dim(UnoiseVOpen_order)

DM["Dada_Bact194", "Deblur_Bact194"]
DavDe_ordered["Dada_Bact194", "Deblur_Bact194"]
DM["Unoise_Bact1z2z2","Dada_Bact1z2z2"]
UvDa["Unoise_Bact1z2z2",which(rownames(UvDa) == "Unoise_Bact1z2z2")]

DM["Unoise_Bact1z2z2","Deblur_Bact1z2z2"]
UvDe_ordered["Unoise_Bact1z2z2",which(rownames(UvDe_ordered) == "Unoise_Bact1z2z2")]

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))

Unifrac_boxplot <- boxplot(diag(DavDe_ordered),
                           diag(UvDa),
                           diag(UvDe_ordered),
                           diag(DadaVOpenorder),
                           diag(DeblurVOpen_order),
                           diag(UnoiseVOpen_order),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open",  "Deblur_Open", "UNOISE3_Open"),
                           ylab="Weighted UniFrac Distance",
                           cex.axis=0.66,
                           outline=FALSE)

Blueberry_box <- recordPlot(Unifrac_boxplot)

#ordination
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

unifrac_combined_blueberry <- sample_cord[,c("V2", "V3", "V4")]
colnames(unifrac_combined_blueberry) <- c("PC1", "PC2", "PC3")
rownames(unifrac_combined_blueberry) <- gsub("001101", "", rownames(unifrac_combined_blueberry))
unifrac_combined_blueberry$sample <- as.factor(gsub("^.+_", "", rownames(unifrac_combined_blueberry)))
unifrac_combined_blueberry$Pipeline <- factor(gsub("_.+$", "", rownames(unifrac_combined_blueberry)))
unifrac_combined_blueberry$soil <- NA
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% rhizosphere_samples), "soil"] <- "Rhizosphere"
#unifrac_combined_blueberry[which(unifrac_combined_blueberry$sample %in% bulk_samples), "soil"] <- "Bulk"


#remove noPos samples
unifrac_combined_blueberry_temp <- unifrac_combined_blueberry[-which(unifrac_combined_blueberry$Pipeline %in% "Deblur-noPos"),]
#remove samples not included in the others
samples_to_keep <- unifrac_combined_blueberry_temp[unifrac_combined_blueberry_temp$Pipeline == "DADA2", "sample"]
samples_to_keep <- as.character(levels(samples_to_keep))


unifrac_combined_filtered <- unifrac_combined_blueberry_temp[unifrac_combined_blueberry_temp$sample %in% samples_to_keep,]
unifrac_combined_filtered$Pipeline <- gsub("Dada", "DADA2", unifrac_combined_filtered$Pipeline)
unifrac_combined_filtered$Pipeline <- gsub("open-ref", "OTU", unifrac_combined_filtered$Pipeline)
unifrac_combined_filtered$Pipeline <- gsub("Unoise", "UNOISE3", unifrac_combined_filtered$Pipeline)

unifrac_plot <- ggplot(data=unifrac_combined_filtered, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(unifrac_combined_filtered$sample)) {
  unifrac_plot <- unifrac_plot + geom_line(data = unifrac_combined_filtered[which(unifrac_combined_filtered$sample == sample),], aes(x = PC1, y = PC2))
}

plot(unifrac_plot)



################ BISCUIT data now
setwd("/home/jacob/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/fixed_combined/pplacer_distances/")

#read in proportions for the PC
Bis_proportions <- read.table("weighted_PC_cords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
Bis_SampleNum <- read.table("weighted_PC_cords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
Bis_sample_cord <- read.table("weighted_PC_cords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = Bis_SampleNum[[2]],
                          header=FALSE, row.names=1)

Bis_unifrac_combined_biscuit <- Bis_sample_cord[,c("V2", "V3", "V4")]
colnames(Bis_unifrac_combined_biscuit) <- c("PC1", "PC2", "PC3")
Bis_unifrac_combined_biscuit$sample <- as.factor(gsub("^.+_", "", rownames(Bis_unifrac_combined_biscuit)))
Bis_unifrac_combined_biscuit$Pipeline <- factor(gsub("_.+$", "", rownames(Bis_unifrac_combined_biscuit)))
Bis_unifrac_combined_biscuit$soil <- NA
Bis_unifrac_combined_biscuit$Pipeline <- gsub("Dada", "DADA2", Bis_unifrac_combined_biscuit$Pipeline)
Bis_unifrac_combined_biscuit$Pipeline <- gsub("Unoise", "UNOISE3", Bis_unifrac_combined_biscuit$Pipeline)
Bis_unifrac_combined_biscuit$Pipeline <- gsub("open-ref", "OTU", Bis_unifrac_combined_biscuit$Pipeline)
Bis_unifrac_combined_biscuit$sample <- gsub("001101", "", Bis_unifrac_combined_biscuit$sample)
rownames(Bis_unifrac_combined_biscuit) <- gsub("001101", "", rownames(Bis_unifrac_combined_biscuit))

#filter out samples that were removed by rarfaction 

samples_to_keep <- Bis_unifrac_combined_biscuit[Bis_unifrac_combined_biscuit$Pipeline == "DADA2", "sample"]

Bis_unifrac_filtered <- Bis_unifrac_combined_biscuit[Bis_unifrac_combined_biscuit$sample %in% samples_to_keep,]
Bis_unifrac_filtered <- Bis_unifrac_filtered[-which(Bis_unifrac_filtered$Pipeline=="Deblur-noPos"),]
Bis_unifrac_filtered$sample <- as.factor(Bis_unifrac_filtered$sample)

Bis_unifrac_plot <- ggplot(data=Bis_unifrac_filtered, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(Bis_proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(Bis_proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)


for (sample in levels(Bis_unifrac_filtered$sample)) {
  Bis_unifrac_plot <- Bis_unifrac_plot + geom_line(data = Bis_unifrac_filtered[which(Bis_unifrac_filtered$sample == sample),], aes(x = PC1, y = PC2))
}

plot(Bis_unifrac_plot)


#Biscuit Boxplot


Bis_DM <- read.table("weighted_unifrac_rare_Combined_biscuit_taxa_fixed.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
colnames(Bis_DM) <- gsub("001101", "", colnames(Bis_DM))
rownames(Bis_DM) <- gsub("001101", "", rownames(Bis_DM))

# Read in blueberry metadata.
#map <- read.table("../../map.txt", sep="\t", header=T, stringsAsFactors = FALSE) 

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
Bis_unoiseVDada <- Bis_DM[grepl("Unoise*", names(Bis_DM)) ,]
Bis_unoiseVDada <- Bis_unoiseVDada[, grepl("Dada*", names(Bis_unoiseVDada))]
Bis_UvDa <- data.matrix(Bis_unoiseVDada)
# Already in order so don't need to re-order!

# Unoise vs Deblur
Bis_unoiseVDeblur <- Bis_DM[grepl("Unoise*", names(Bis_DM)),]
Bis_unoiseVDeblur <- Bis_unoiseVDeblur[, grepl("Deblur*", names(Bis_unoiseVDeblur))]
Bis_UvDe <- data.matrix(Bis_unoiseVDeblur)

Bis_UvDe_ordered <- reorder_pipeline_samples(dm_input=Bis_UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

# Unoise vs Deblur
Bis_dadaVdeblur <- Bis_DM[grepl("Dada*", names(Bis_DM)) ,]
Bis_dadaVdeblur <- Bis_dadaVdeblur[,grepl("Deblur*", names(Bis_dadaVdeblur))]
Bis_DavDe <- data.matrix(Bis_dadaVdeblur)

Bis_DavDe_ordered <- reorder_pipeline_samples(dm_input=Bis_DavDe, row_prefix="Dada_", col_prefix="Deblur_")


# Sanity checks:
dim(Bis_UvDa)
dim(Bis_UvDe_ordered)
dim(Bis_DavDe_ordered)


#filter out samples in Open that passed rarefaction but others did not.
Bis_samples_to_keep <- colnames(Bis_DM)[grep("Deblur_", colnames(Bis_DM))]
Bis_samples_to_keep <- gsub("Deblur_", "", Bis_samples_to_keep)
Bis_DM_filt <- Bis_DM[grep(paste(Bis_samples_to_keep, collapse = "|"), rownames(Bis_DM)), grep(paste(Bis_samples_to_keep, collapse = "|"), colnames(Bis_DM))]

#UnoiseVOpen
Bis_UnoiseVOpen <- Bis_DM_filt[grep("Unoise", names(Bis_DM_filt)), grep("open-ref", names(Bis_DM_filt))]
Bis_UnoiseVOpen_order  <- reorder_pipeline_samples(dm_input = data.matrix(Bis_UnoiseVOpen), row_prefix = "Unoise__", col_prefix = "open-ref_")

#DadaVOpen
Bis_DadaVOpen <- Bis_DM_filt[grep("Dada_", names(Bis_DM_filt)), grep("open-ref", names(Bis_DM_filt))]
Bis_DadaVOpen_order <- reorder_pipeline_samples(dm_input = data.matrix(Bis_DadaVOpen), row_prefix = "Dada_", col_prefix = "open-ref_")

#DeblurVOpen
Bis_DeblurVOpen <- Bis_DM_filt[grep("Deblur_", names(Bis_DM_filt)), grep("open-ref_", names(Bis_DM_filt))]
Bis_DeblurVOpen_order <- reorder_pipeline_samples(dm_input = data.matrix(Bis_DeblurVOpen), row_prefix = "Deblur_", col_prefix = "open-ref_")

#Deblur V Deblur no Pos
Bis_DeblurVnoPos <- Bis_DM_filt[grep("Deblur_", names(Bis_DM_filt)), grep("Deblur-noPos_", names(Bis_DM_filt))]
Bis_DeblurVnoPos_order <- reorder_pipeline_samples(dm_input = data.matrix(Bis_DeblurVnoPos), row_prefix = "Deblur_", col_prefix = 'Deblur-noPos_')

gsub("Dada_", "", colnames(Bis_UvDa)) == gsub("Unoise__", "", rownames(Bis_UvDa))
gsub("Deblur_", "", colnames(Bis_UvDe_ordered)) == gsub("Unoise__", "", rownames(Bis_UvDe_ordered))
gsub("Deblur_", "", colnames(Bis_DavDe_ordered)) == gsub("Dada_", "", rownames(Bis_DavDe_ordered))
gsub("Dada_", "", rownames(Bis_DadaVOpen_order)) == gsub("open-ref_", "", colnames(Bis_DadaVOpen_order))
gsub("Deblur_", "", rownames(Bis_DeblurVOpen_order)) == gsub("open-ref_", "", colnames(Bis_DeblurVOpen_order))
gsub("Unoise__", "", rownames(Bis_UnoiseVOpen_order)) == gsub("open-ref_", "", colnames(Bis_UnoiseVOpen_order))
gsub("Deblur_", "", rownames(Bis_DeblurVnoPos_order)) == gsub("Deblur-noPos_", "", colnames(Bis_DeblurVnoPos_order))


#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Bis_Unifrac_boxplot <- boxplot(diag(Bis_DavDe_ordered),
                           diag(Bis_UvDa),
                           diag(Bis_UvDe_ordered),
                           diag(Bis_DadaVOpen_order),
                           diag(Bis_DeblurVOpen_order),
                           diag(Bis_UnoiseVOpen_order),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_OTU", "Deblur_OTU", "UNOISE3_OTU"),
                           ylab="Weighted UniFrac Distance",
                           cex.axis=0.66,
                           outline=FALSE)

#Used to make figure 4
Bis_Unifrac_box <- recordPlot(Bis_Unifrac_boxplot)





#### Exercise Data

setwd("/home/jacob/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/fixed_combined/pplacer_distances/")

### Plot weighted UniFrac data

#read in proportions for the PC
E_proportions <- read.table("weighted_PC_cords.txt", sep = "\t", nrow=1, skip=4)

#read in the number oof samples
E_SampleNum <- read.table("weighted_PC_cords.txt", sep = "\t", nrow=1)

#Get coordinates for samples.
E_sample_cord <- read.table("weighted_PC_cords.txt", 
                          sep ="\t", skip = 9, 
                          nrow = E_SampleNum[[2]],
                          header=FALSE, row.names=1)

E_unifrac_combined_Exercise <- E_sample_cord[,c("V2", "V3", "V4")]
colnames(E_unifrac_combined_Exercise) <- c("PC1", "PC2", "PC3")
E_unifrac_combined_Exercise$sample <- as.factor(gsub("^.+_", "", rownames(E_unifrac_combined_Exercise)))
E_unifrac_combined_Exercise$Pipeline <- factor(gsub("_.+$", "", rownames(E_unifrac_combined_Exercise)))
E_unifrac_combined_Exercise$Pipeline <- gsub("Dada", "DADA2", E_unifrac_combined_Exercise$Pipeline)
E_unifrac_combined_Exercise$Pipeline <- gsub("Unoise", "UNOISE3", E_unifrac_combined_Exercise$Pipeline)
E_unifrac_combined_Exercise$Pipeline <- gsub("open-ref", "OTU", E_unifrac_combined_Exercise$Pipeline)
E_unifrac_combined_rm_noPos <- E_unifrac_combined_Exercise[-which(E_unifrac_combined_Exercise$Pipeline=="Deblur-noPos"),]

E_samples_to_keep <- E_unifrac_combined_rm_noPos[E_unifrac_combined_rm_noPos$Pipeline=="Deblur", "sample"]

E_unifrac_combined_filtered <- E_unifrac_combined_rm_noPos[which(E_unifrac_combined_rm_noPos$sample %in% E_samples_to_keep),]
#need to filter the sample list so that there ar enot single samples... 

E_unifrac_plot <- ggplot(data=E_unifrac_combined_filtered, aes(PC1, PC2)) +
  geom_point(aes(fill=sample, size=1.5, shape=Pipeline)) +
  theme_minimal() +
  xlab(paste("PC1 (",round(E_proportions$V1, 2)*100,"%)", sep="")) +
  ylab(paste("PC2 (",round(E_proportions$V2, 2)*100,"%)", sep="")) +
  scale_shape_manual(values = c(21, 22, 23, 24 ,25)) +
  guides(size=FALSE, fill=FALSE) +
  scale_fill_manual(values=diff_col)

for (sample in levels(E_unifrac_combined_filtered$sample)) {
  E_unifrac_plot <- E_unifrac_plot + geom_line(data = E_unifrac_combined_filtered[which(E_unifrac_combined_filtered$sample == sample),], aes(x = PC1, y = PC2))
}

plot(E_unifrac_plot)




### Exercise boxplot

E_DM <- read.table("weighted_unifrac_rare_Combined_Exercise_fixed.txt", sep="\t", header=TRUE, row.names=1)



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
E_unoiseVDada <- E_DM[grepl("Unoise*", names(E_DM)) ,]
E_unoiseVDada <- E_unoiseVDada[, grepl("Dada*", names(E_unoiseVDada))]
E_UvDa <- data.matrix(E_unoiseVDada)
# Already in order so don't need to re-order!

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


#sample to keep
E_samples_to_keep <- colnames(E_DM)[grep("Deblur_", colnames(E_DM))]
E_samples_to_keep <- gsub("Deblur_", "open-ref_", E_samples_to_keep)

#Vs Open figures
E_DadaVOpen <- E_DM[grep("open-ref_", rownames(E_DM)), grep("Dada_", names(E_DM))]
#filter out samples
E_DadaVOpen <- E_DadaVOpen[rownames(E_DadaVOpen) %in% E_samples_to_keep, ]
E_DadaVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(E_DadaVOpen), row_prefix = "open-ref_", col_prefix = "Dada_")

E_DeblurVOpen <- E_DM[grep("open-ref_", rownames(E_DM)), grep("Deblur_", names(E_DM))]
#filter out samples
E_DeblurVOpen <- E_DeblurVOpen[rownames(E_DeblurVOpen) %in% E_samples_to_keep,]
E_DeblurVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(E_DeblurVOpen), row_prefix = "open-ref_", col_prefix = "Deblur_")

E_UnoiseVOpen <- E_DM[grep("open-ref_", rownames(E_DM)), grep("Unoise_", names(E_DM))]
E_UnoiseVOpen <- E_UnoiseVOpen[rownames(E_UnoiseVOpen) %in% E_samples_to_keep, ]
E_UnoiseVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(E_UnoiseVOpen), row_prefix = "open-ref_", col_prefix = "Unoise__")

# Sanity checks:
dim(E_UvDa)
dim(E_UvDe_ordered)
dim(E_DavDe_ordered)
dim(E_DadaVOpen_sorted)
dim(E_DeblurVOpen_sorted)
dim(E_UnoiseVOpen_sorted)


#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
E_Unifrac_boxplot <- boxplot(diag(E_DavDe_ordered),
                           diag(E_UvDa),
                           diag(E_UvDe_ordered),
                           diag(E_DadaVOpen_sorted),
                           diag(E_DeblurVOpen_sorted),
                           diag(E_UnoiseVOpen_sorted),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_OTU", "Deblur_OTU", "UNOISE3_OTU"),
                           ylab="Weighted UniFrac Distance",
                           cex.axis=0.66,
                           outline=FALSE)

#Used to make figure 4
E_Unifrac_box <- recordPlot(E_Unifrac_boxplot)





library(cowplot)


weighted_plot_fig2 <- plot_grid(Blueberry_box, Bis_Unifrac_box, E_Unifrac_box, unifrac_plot, Bis_unifrac_plot, E_unifrac_plot, labels="AUTO", rel_heights = c(1,2))
weighted_plot_fig2




