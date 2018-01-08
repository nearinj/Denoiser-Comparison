#Used to make the Weighted Unifrac Box plots for the blueberry data

library(ggplot2)

setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/plots/bdiv-out/")
par=(mfrow=c(2,2))

# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
DM <- read.table("weighted_unifrac_dm.txt", sep="\t", header=TRUE, row.names=1)

# Read in blueberry metadata.
bluberry_map <- read.table("../../map_blueberry_merged.csv", sep="\t", header=T, stringsAsFactors = FALSE,) 

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
        names=c("Dada2_Deblur", "Dada2_Unoise3", "Unoise3_Deblur"),
        ylab="Weighted UniFrac Distance", 
        outline=FALSE)

#Used to make figure 4
Unifrac_box <- recordPlot(Unifrac_boxplot)

Rhizosphere <- bluberry_map[grepl("Rhizosphere", bluberry_map$Description) ,]
Bulk <- bluberry_map[grepl("Bulk", bluberry_map$Description),]



dadaVdadaBulk <- dadaVdada[ rownames(dadaVdada) %in% Bulk$SampleID, colnames(dadaVdada) %in% Bulk$SampleID]
dadaVdadaRhizo <- dadaVdada[ rownames(dadaVdada) %in% Rhizosphere$SampleID, colnames(dadaVdada) %in% Rhizosphere$SampleID]

identical(colnames(dadaVdadaRhizo), rownames(dadaVdadaRhizo))
identical(colnames(dadaVdadaBulk), rownames(dadaVdadaBulk))

deblurVdeblurBulk <- deblurVdeblur[ rownames(deblurVdeblur) %in% Bulk$SampleID, colnames(deblurVdeblur) %in% Bulk$SampleID]
deblurVdeblurRhizo <- deblurVdeblur[ rownames(deblurVdeblur) %in% Rhizosphere$SampleID, colnames(deblurVdeblur) %in% Rhizosphere$SampleID]

identical(colnames(deblurVdeblurRhizo), rownames(deblurVdeblurRhizo))
identical(colnames(deblurVdeblurBulk), rownames(deblurVdeblurBulk))


unoiseVunoiseBulk <- unoiseVunoise[ rownames(unoiseVunoise) %in% Bulk$SampleID, colnames(unoiseVunoise) %in% Bulk$SampleID]
unoiseVunoiseRhizo <- unoiseVunoise[ rownames(unoiseVunoise) %in% Rhizosphere$SampleID, colnames(unoiseVunoise) %in% Rhizosphere$SampleID]

identical(colnames(unoiseVunoiseBulk), rownames(unoiseVunoiseBulk))
identical(colnames(unoiseVunoiseRhizo), rownames(unoiseVunoiseRhizo))

#box plot of intra sample distances broken by rhizosphere and bulk!
boxplot(diag(DavDe_ordered),
        diag(UvDa),
        diag(UvDe_ordered),
        deblurVdeblurBulk[upper.tri(deblurVdeblurBulk, diag=FALSE)],
        dadaVdadaBulk[upper.tri(dadaVdadaBulk, diag=FALSE)],
        unoiseVunoiseBulk[upper.tri(unoiseVunoiseBulk, diag=FALSE)],
        deblurVdeblurRhizo[upper.tri(deblurVdeblurRhizo, diag=FALSE)],
        dadaVdadaRhizo[upper.tri(dadaVdadaRhizo, diag=FALSE)],
        unoiseVunoiseRhizo[upper.tri(unoiseVunoiseRhizo, diag=FALSE)],
        names=c("Dada_Deblur", "Unoise_Dada","Unoise_Deblur","Deblur_Bulk", "Dada_Bulk", "Unoise_Bulk", "Deblur_Rhizo", "Dada_Rhizo", "Unoise_Rhizo"),
        outline = FALSE,
        ylab="Weighted UniFrac Distance"
  )

