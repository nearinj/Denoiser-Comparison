#Used to make the Weighted Unifrac Box plots for the blueberry data

library(ggplot2)

setwd("/home/jacob/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/fixed_combined/pplacer_distances/")
par=(mfrow=c(2,2))

# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
DM <- read.table("weighted_unifrac_rare_Combined_Exercise_fixed.txt", sep="\t", header=TRUE, row.names=1)



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

UvDe_ordered <- reorder_pipeline_samples(dm_input=UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

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

#sample to keep
samples_to_keep <- colnames(DM)[grep("Deblur_", colnames(DM))]
samples_to_keep <- gsub("Deblur_", "open-ref_", samples_to_keep)

#Vs Open figures
DadaVOpen <- DM[grep("open-ref_", rownames(DM)), grep("Dada_", names(DM))]
#filter out samples
DadaVOpen <- DadaVOpen[rownames(DadaVOpen) %in% samples_to_keep, ]
DadaVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(DadaVOpen), row_prefix = "open-ref_", col_prefix = "Dada_")

DeblurVOpen <- DM[grep("open-ref_", rownames(DM)), grep("Deblur_", names(DM))]
#filter out samples
DeblurVOpen <- DeblurVOpen[rownames(DeblurVOpen) %in% samples_to_keep,]
DeblurVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(DeblurVOpen), row_prefix = "open-ref_", col_prefix = "Deblur_")

UnoiseVOpen <- DM[grep("open-ref_", rownames(DM)), grep("Unoise_", names(DM))]
UnoiseVOpen <- UnoiseVOpen[rownames(UnoiseVOpen) %in% samples_to_keep, ]
UnoiseVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(UnoiseVOpen), row_prefix = "open-ref_", col_prefix = "Unoise__")

# Sanity checks:
dim(UvDa)
dim(UvDe_ordered)
dim(DavDe_ordered)
dim(DadaVOpen_sorted)
dim(DeblurVOpen_sorted)
dim(UnoiseVOpen_sorted)


#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(DavDe_ordered),
        diag(UvDa),
        diag(UvDe_ordered),
        diag(DadaVOpen_sorted),
        diag(DeblurVOpen_sorted),
        diag(UnoiseVOpen_sorted),
        names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open", "Deblur_Open", "UNOISE3_open"),
        ylab="Weighted UniFrac Distance",
        cex.axis=0.66,
        outline=FALSE)

#Used to make figure 4
Unifrac_box <- recordPlot(Unifrac_boxplot)

saveRDS(Unifrac_box, "../weighted_unifrac_boxplot.rds")

#################################################################unweighted unifrac##############################################################

U_DM <- read.table("unweighted_unifrac_rare_Combined_Exercise_fixed.txt", sep="\t", header=TRUE, row.names=1)

# Unoise vs Dada
U_unoiseVDada <- U_DM[grepl("Unoise*", names(U_DM)) ,]
U_unoiseVDada <- U_unoiseVDada[, grepl("Dada*", names(U_unoiseVDada))]
U_UvDa <- data.matrix(U_unoiseVDada)
# Already in order so don't need to re-order!

# Unoise vs Deblur
U_unoiseVDeblur <- U_DM[grepl("Unoise*", names(U_DM)),]
U_unoiseVDeblur <- U_unoiseVDeblur[, grepl("Deblur*", names(U_unoiseVDeblur))]
U_UvDe <- data.matrix(U_unoiseVDeblur)

U_UvDe_ordered <- reorder_pipeline_samples(dm_input=U_UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

# Unoise vs Deblur
U_dadaVdeblur <- U_DM[grepl("Dada*", names(U_DM)) ,]
U_dadaVdeblur <- U_dadaVdeblur[,grepl("Deblur*", names(U_dadaVdeblur))]
U_DavDe <- data.matrix(U_dadaVdeblur)

U_DavDe_ordered <- reorder_pipeline_samples(dm_input=U_DavDe, row_prefix="Dada_", col_prefix="Deblur_")

#sample to keep
U_samples_to_keep <- colnames(U_DM)[grep("Deblur_", colnames(U_DM))]
U_samples_to_keep <- gsub("Deblur_", "open-ref_", U_samples_to_keep)

#Vs Open figures
U_DadaVOpen <- U_DM[grep("open-ref_", rownames(U_DM)), grep("Dada_", names(U_DM))]
#filter out samples
U_DadaVOpen <- U_DadaVOpen[rownames(U_DadaVOpen) %in% U_samples_to_keep, ]
U_DadaVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(U_DadaVOpen), row_prefix = "open-ref_", col_prefix = "Dada_")

U_DeblurVOpen <- U_DM[grep("open-ref_", rownames(U_DM)), grep("Deblur_", names(U_DM))]
#filter out samples
U_DeblurVOpen <- U_DeblurVOpen[rownames(U_DeblurVOpen) %in% U_samples_to_keep,]
U_DeblurVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(U_DeblurVOpen), row_prefix = "open-ref_", col_prefix = "Deblur_")

U_UnoiseVOpen <- U_DM[grep("open-ref_", rownames(U_DM)), grep("Unoise_", names(U_DM))]
U_UnoiseVOpen <- U_UnoiseVOpen[rownames(U_UnoiseVOpen) %in% U_samples_to_keep, ]
U_UnoiseVOpen_sorted <- reorder_pipeline_samples(dm_input = as.matrix(U_UnoiseVOpen), row_prefix = "open-ref_", col_prefix = "Unoise__")

# Sanity checks:
dim(U_UvDa)
dim(U_UvDe_ordered)
dim(U_DavDe_ordered)
dim(U_DadaVOpen_sorted)
dim(U_DeblurVOpen_sorted)
dim(U_UnoiseVOpen_sorted)


#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
U_Unifrac_boxplot <- boxplot(diag(U_DavDe_ordered),
                           diag(U_UvDa),
                           diag(U_UvDe_ordered),
                           diag(U_DadaVOpen_sorted),
                           diag(U_DeblurVOpen_sorted),
                           diag(U_UnoiseVOpen_sorted),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open", "Deblur_Open", "UNOISE3_open"),
                           ylab="Unweighted UniFrac Distance",
                           cex.axis=0.66,
                           outline=FALSE)

#Used to make figure 4
U_Unifrac_box <- recordPlot(U_Unifrac_boxplot)
Exercise_U_Unifrac_box <- U_Unifrac_box
