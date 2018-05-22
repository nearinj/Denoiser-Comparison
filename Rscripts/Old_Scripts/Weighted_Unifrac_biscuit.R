#Used to make the Weighted Unifrac Box plots for the blueberry data

library(ggplot2)

setwd("/home/jacob/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/fixed_combined/pplacer_distances/")
par=(mfrow=c(2,2))

# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
DM <- read.table("weighted_unifrac_rare_Combined_biscuit_taxa_fixed.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
colnames(DM) <- gsub("001101", "", colnames(DM))
rownames(DM) <- gsub("001101", "", rownames(DM))

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
deblurVdeblur <- DM[grep("Deblur", names(DM)) ,grep("Deblur", names(DM))]
DevDe <- data.matrix(deblurVdeblur)
DevDe_ordered <- reorder_pipeline_samples(dm_input=DevDe, row_prefix="Deblur_", col_prefix="Deblur_")

# unoise vs unoise
unoiseVunoise <- DM[grep("Unoise", names(DM)) ,grep("Unoise", names(DM))]
UvU <- data.matrix(unoiseVunoise)
UvU_ordered <- reorder_pipeline_samples(dm_input=UvU, row_prefix="Unoise_", col_prefix="Unoise_")


# Sanity checks:
dim(UvDa)
dim(UvDe_ordered)
dim(DavDe_ordered)


#filter out samples in Open that passed rarefaction but others did not.
samples_to_keep <- colnames(DM)[grep("Deblur_", colnames(DM))]
samples_to_keep <- gsub("Deblur_", "", samples_to_keep)
DM_filt <- DM[grep(paste(samples_to_keep, collapse = "|"), rownames(DM)), grep(paste(samples_to_keep, collapse = "|"), colnames(DM))]

#UnoiseVOpen
UnoiseVOpen <- DM_filt[grep("Unoise", names(DM_filt)), grep("open-ref", names(DM_filt))]
UnoiseVOpen_order  <- reorder_pipeline_samples(dm_input = data.matrix(UnoiseVOpen), row_prefix = "Unoise__", col_prefix = "open-ref_")

#DadaVOpen
DadaVOpen <- DM_filt[grep("Dada_", names(DM_filt)), grep("open-ref", names(DM_filt))]
DadaVOpen_order <- reorder_pipeline_samples(dm_input = data.matrix(DadaVOpen), row_prefix = "Dada_", col_prefix = "open-ref_")

#DeblurVOpen
DeblurVOpen <- DM_filt[grep("Deblur_", names(DM_filt)), grep("open-ref_", names(DM_filt))]
DeblurVOpen_order <- reorder_pipeline_samples(dm_input = data.matrix(DeblurVOpen), row_prefix = "Deblur_", col_prefix = "open-ref_")

#Deblur V Deblur no Pos
DeblurVnoPos <- DM_filt[grep("Deblur_", names(DM_filt)), grep("Deblur-noPos_", names(DM_filt))]
DeblurVnoPos_order <- reorder_pipeline_samples(dm_input = data.matrix(DeblurVnoPos), row_prefix = "Deblur_", col_prefix = 'Deblur-noPos_')

gsub("Dada_", "", colnames(UvDa)) == gsub("Unoise__", "", rownames(UvDa))
gsub("Deblur_", "", colnames(UvDe_ordered)) == gsub("Unoise__", "", rownames(UvDe_ordered))
gsub("Deblur_", "", colnames(DavDe_ordered)) == gsub("Dada_", "", rownames(DavDe_ordered))
gsub("Dada_", "", rownames(DadaVOpen_order)) == gsub("open-ref_", "", colnames(DadaVOpen_order))
gsub("Deblur_", "", rownames(DeblurVOpen_order)) == gsub("open-ref_", "", colnames(DeblurVOpen_order))
gsub("Unoise__", "", rownames(UnoiseVOpen_order)) == gsub("open-ref_", "", colnames(UnoiseVOpen_order))
gsub("Deblur_", "", rownames(DeblurVnoPos_order)) == gsub("Deblur-noPos_", "", colnames(DeblurVnoPos_order))


#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(DavDe_ordered),
        diag(UvDa),
        diag(UvDe_ordered),
        diag(DadaVOpen_order),
        diag(DeblurVOpen_order),
        diag(UnoiseVOpen_order),
        names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open", "Deblur_Open", "UNOISE3_Open"),
        ylab="Weighted UniFrac Distance",
        cex.axis=0.66,
        outline=FALSE)

#Used to make figure 4
Unifrac_box <- recordPlot(Unifrac_boxplot)

################################################################unweighted unifrac##############################################

#get unweighted distances between samples
# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
U_DM <- read.table("unweighted_unifrac_rare_Combined_biscuit_taxa_fixed.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
colnames(U_DM) <- gsub("001101", "", colnames(U_DM))
rownames(U_DM) <- gsub("001101", "", rownames(U_DM))

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



# Sanity checks:
dim(UvDa)
dim(UvDe_ordered)
dim(DavDe_ordered)


#filter out samples in Open that passed rarefaction but others did not.
U_samples_to_keep <- colnames(U_DM)[grep("Deblur_", colnames(U_DM))]
U_samples_to_keep <- gsub("Deblur_", "", U_samples_to_keep)
U_DM_filt <- U_DM[grep(paste(U_samples_to_keep, collapse = "|"), rownames(U_DM)), grep(paste(U_samples_to_keep, collapse = "|"), colnames(U_DM))]

#UnoiseVOpen
U_UnoiseVOpen <- U_DM_filt[grep("Unoise", names(U_DM_filt)), grep("open-ref", names(U_DM_filt))]
U_UnoiseVOpen_order  <- reorder_pipeline_samples(dm_input = data.matrix(U_UnoiseVOpen), row_prefix = "Unoise__", col_prefix = "open-ref_")

#DadaVOpen
U_DadaVOpen <- U_DM_filt[grep("Dada_", names(U_DM_filt)), grep("open-ref", names(U_DM_filt))]
U_DadaVOpen_order <- reorder_pipeline_samples(dm_input = data.matrix(U_DadaVOpen), row_prefix = "Dada_", col_prefix = "open-ref_")

#DeblurVOpen
U_DeblurVOpen <- U_DM_filt[grep("Deblur_", names(U_DM_filt)), grep("open-ref_", names(U_DM_filt))]
U_DeblurVOpen_order <- reorder_pipeline_samples(dm_input = data.matrix(U_DeblurVOpen), row_prefix = "Deblur_", col_prefix = "open-ref_")


gsub("Dada_", "", colnames(U_UvDa)) == gsub("Unoise__", "", rownames(U_UvDa))
gsub("Deblur_", "", colnames(U_UvDe_ordered)) == gsub("Unoise__", "", rownames(U_UvDe_ordered))
gsub("Deblur_", "", colnames(U_DavDe_ordered)) == gsub("Dada_", "", rownames(U_DavDe_ordered))
gsub("Dada_", "", rownames(U_DadaVOpen_order)) == gsub("open-ref_", "", colnames(U_DadaVOpen_order))
gsub("Deblur_", "", rownames(U_DeblurVOpen_order)) == gsub("open-ref_", "", colnames(U_DeblurVOpen_order))
gsub("Unoise__", "", rownames(U_UnoiseVOpen_order)) == gsub("open-ref_", "", colnames(U_UnoiseVOpen_order))



#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(U_DavDe_ordered),
                           diag(U_UvDa),
                           diag(U_UvDe_ordered),
                           diag(U_DadaVOpen_order),
                           diag(U_DeblurVOpen_order),
                           diag(U_UnoiseVOpen_order),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open", "Deblur_Open", "UNOISE3_Open"),
                           ylab="Weighted UniFrac Distance",
                           cex.axis=0.66,
                           outline=FALSE)

#Used to make figure 4
biscuit_unweighted_Unifrac_box <- recordPlot(Unifrac_boxplot)

