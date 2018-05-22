#Used to make the Bray_Curtis Box plots for the blueberry data

library(ggplot2)
setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/final_combined/taxa_summary/bray_curtis/")

par=(mfrow=c(2,2))
# Read in bray-curtis distance matrix (at genus level) of samples from all 3 pipelines
DM <- read.table("bray_curtis_rare_CombinedV2_Blueberry_taxa_L6.txt", sep="\t", header=TRUE, row.names=1, check.names = F)

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
unoiseVDada <- DM[grepl("Unoise*", names(DM)) ,]
unoiseVDada <- unoiseVDada[, grepl("Dada*", names(unoiseVDada))]
UvDa <- data.matrix(unoiseVDada)
# Already in order so don't need to re-order this matrix.

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
unoiseVunoise <- unoiseVunoise[,grep("Unoise", names(unoiseVunoise))]
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

#DeblurVnoPos
DeblurVnoPos <- as.matrix(DM[grep("Deblur_", names(DM)), grep("Deblur-noPos_", names(DM))])
colnames(DeblurVnoPos) <- gsub("001101", "", colnames(DeblurVnoPos))
DeVnoPos <- reorder_pipeline_samples(DeblurVnoPos, row_prefix = "Deblur_", col_prefix = "Deblur-noPos_")

#DeblurVOpen
DeblurVOpen <- as.matrix(DM[grep("Deblur_", names(DM)), grep("open-ref_", names(DM))])
DeVOP <- reorder_pipeline_samples(DeblurVOpen, "Deblur_", "open-ref_")

#DadaVOpen
DadaVOpen <- as.matrix(DM[grep("Dada_", names(DM)), grep("open-ref_", names(DM))])
DaVOP <- reorder_pipeline_samples(DadaVOpen, "Dada_", "open-ref_")

#UnoiseVOpen
UnoiseVOpen <- as.matrix(DM[grep("Unoise_", names(DM)), grep("open-ref_", names(DM))])
UVOP <- reorder_pipeline_samples(UnoiseVOpen, "Unoise_", "open-ref_")

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Bray_Curits_Boxplot <- boxplot(diag(DavDe_ordered),
        diag(UvDa),
        diag(UvDe_ordered),
        diag(DaVOP),
        diag(DeVOP),
        diag(UVOP),
        names=c("Dada2_Deblur", "Dada2_Unoise3", "Unoise3_Deblur","Open-Ref_DADA2", "Open-Ref_Deblur", "Open-Ref_UNOISE3"),
        ylab="Bray-Curtis Distance",
        cex.axis=0.5,
        outline=FALSE)


#Used to make figure 4
Bray_Curtis_box <- recordPlot(Bray_Curits_Boxplot)










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
        ylab="Bray Curtis Distance"
  )

