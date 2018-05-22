#Used to make the Weighted Unifrac Box plots for the blueberry data

library(ggplot2)

setwd("/home/jacob/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/plots/bdiv-new/")
par=(mfrow=c(2,2))

# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
DM <- read.table("weighted_unifrac_dm.txt", sep="\t", header=TRUE, row.names=1)

# Read in blueberry metadata.
map <- read.table("../../map.txt", sep="\t", header=T, stringsAsFactors = FALSE,) 

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





#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(DavDe_ordered),
        diag(UvDa),
        diag(UvDe_ordered), 
        names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur"),
        ylab="Weighted UniFrac Distance", 
        outline=FALSE)

#Used to make figure 4
Unifrac_box <- recordPlot(Unifrac_boxplot)



#################################################################unweighted unifrac##############################################################

# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
unweighted_DM <- read.table("../bdiv-new/unweighted_unifrac_dm.txt", sep="\t", header=TRUE, row.names=1)


# Unoise vs Dada
unweighted_unoiseVDada <- unweighted_DM[grepl("Unoise*", names(unweighted_DM)) ,]
unweighted_unoiseVDada <- unweighted_unoiseVDada[, grepl("Dada*", names(unweighted_unoiseVDada))]
unweighted_UvDa <- data.matrix(unweighted_unoiseVDada)
# Already in order so don't need to re-order!

# Unoise vs Deblur
unweighted_unoiseVDeblur <- unweighted_DM[grepl("Unoise*", names(unweighted_DM)),]
unweighted_unoiseVDeblur <- unweighted_unoiseVDeblur[, grepl("Deblur*", names(unweighted_unoiseVDeblur))]
unweighted_UvDe <- data.matrix(unweighted_unoiseVDeblur)

unweighted_UvDe_ordered <- reorder_pipeline_samples(dm_input=unweighted_UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

# Unoise vs Deblur
unweighted_dadaVdeblur <- unweighted_DM[grepl("Dada*", names(unweighted_DM)) ,]
unweighted_dadaVdeblur <- unweighted_dadaVdeblur[,grepl("Deblur*", names(unweighted_dadaVdeblur))]
unweighted_DavDe <- data.matrix(unweighted_dadaVdeblur)

unweighted_DavDe_ordered <- reorder_pipeline_samples(dm_input=unweighted_DavDe, row_prefix="Dada_", col_prefix="Deblur_")

# Dada vs Dada
unweighted_dadaVdada <- unweighted_DM[grep("Dada", names(unweighted_DM)) ,grep("Dada", names(unweighted_DM))]
unweighted_DavDa <- data.matrix(unweighted_dadaVdada)
unweighted_DavDa_ordered <- reorder_pipeline_samples(dm_input=unweighted_DavDa, row_prefix="Dada_", col_prefix="Dada_")

# Deblur vs Deblur
unweighted_deblurVdeblur <- unweighted_DM[grep("Deblur", names(unweighted_DM)) ,grep("Deblur", names(unweighted_DM))]
unweighted_DevDe <- data.matrix(unweighted_deblurVdeblur)
unweighted_DevDe_ordered <- reorder_pipeline_samples(dm_input=unweighted_DevDe, row_prefix="Deblur_", col_prefix="Deblur_")

# unoise vs unoise
unweighted_unoiseVunoise <- unweighted_DM[grep("Unoise", names(unweighted_DM)) ,grep("Unoise", names(unweighted_DM))]
unweighted_UvU <- data.matrix(unweighted_unoiseVunoise)
unweighted_UvU_ordered <- reorder_pipeline_samples(dm_input=unweighted_UvU, row_prefix="Unoise_", col_prefix="Unoise_")


# Sanity checks:
dim(unweighted_UvDa)
dim(unweighted_UvDe_ordered)
dim(unweighted_DavDe_ordered)





#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
unweighted_unifrac_boxplot <- boxplot(diag(unweighted_DavDe_ordered),
                           diag(unweighted_UvDa),
                           diag(unweighted_UvDe_ordered), 
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur"),
                           ylab="Weighted UniFrac Distance", 
                           outline=FALSE)

#make these plots now and then ill be done with the unweighted_unifrac analysis except for the otu picking analysis.....