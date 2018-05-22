#Used to make the Bray_Curtis Box plots for the blueberry data

library(ggplot2)
setwd("/home/jacob/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/final_combined/bray_curtis/distances/")

par=(mfrow=c(2,2))
# Read in bray-curtis distance matrix (at genus level) of samples from all 3 pipelines
DM <- read.table("bray_curtis_rare_CombinedV2_Exercise_taxa_L6.txt", sep="\t", header=TRUE, row.names=1, check.names=F)

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
unoiseVDada <- DM[grepl("Unoise*", names(DM)) ,]

unoiseVDada <- unoiseVDada[, grepl("Dada*", names(unoiseVDada))]
UvDa <- data.matrix(unoiseVDada)
# Already in order so don't need to re-order this matrix.

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
unoiseVunoise <- unoiseVunoise[,grep("Unoise", names(unoiseVunoise))]
UvU <- data.matrix(unoiseVunoise)
UvU_ordered <- reorder_pipeline_samples(dm_input=UvU, row_prefix="Unoise_", col_prefix="Unoise_")


# Sanity checks:
dim(UvDa)
dim(UvDe_ordered)
dim(DavDe_ordered)

#DeblurVnoPos
DeblurVnoPos <- as.matrix(DM[grep("Deblur_", names(DM)), grep("Deblur-noPos_", names(DM))])
colnames(DeblurVnoPos) <- gsub("001101", "", colnames(DeblurVnoPos))
DeVnoPos <- reorder_pipeline_samples(DeblurVnoPos, "Deblur_", "Deblur-noPos_")

#Deblur V Open
DeblurVOpen <- as.matrix(DM[grep("Deblur_", names(DM)), grep("open-ref_", names(DM))])
DeVOP <- reorder_pipeline_samples(DeblurVOpen, "Deblur_", "open-ref_")

#Dada V Open
DadaVOpen <- as.matrix(DM[grep("Dada_", names(DM)), grep("open-ref_", names(DM))])
DaVOP <- reorder_pipeline_samples(DadaVOpen, "Dada_", "open-ref_")

#Unoise V Open

UnoiseVOpen <- as.matrix(DM[grep("Unoise_", names(DM)), grep("open-ref_", names(DM))])
UVOP <- reorder_pipeline_samples(UnoiseVOpen, "Unoise__", "open-ref_")

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Bray_Curits_Boxplot <- boxplot(diag(DavDe_ordered),
        diag(UvDa),
        diag(UvDe_ordered),
        diag(DaVOP),
        diag(DeVOP),
        diag(UVOP),
        names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open", "Debur_Open", "UNOISE3_Open"),
        ylab="Bray-Curtis Distance",
        cex.axis=0.66,
        outline=FALSE)

#Used to make figure 4
Bray_Curtis_box <- recordPlot(Bray_Curits_Boxplot)

