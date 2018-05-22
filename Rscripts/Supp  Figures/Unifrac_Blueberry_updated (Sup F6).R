#Used to make the Weighted Unifrac Box plots for the blueberry data

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


#noPos
rownames(DM) <- gsub("001101", "", rownames(DM))
colnames(DM) <- gsub("001101", "", colnames(DM))


#DadaVnoPOs 
DadaVnoPos <- as.matrix(DM[grep("Dada_", rownames(DM)), grep("Deblur.noPos_", colnames(DM))])
DadaVnoPos_order <- reorder_pipeline_samples(DadaVnoPos, row_prefix = "Dada_", col_prefix = "Deblur.noPos_")

#UnoiseVnoPos
UnoiseVnoPos <- as.matrix(DM[grep("Unoise_", rownames(DM)), grep("Deblur.noPos_", colnames(DM))])
UnoiseVnoPos_order <- reorder_pipeline_samples(UnoiseVnoPos, row_prefix = "Unoise_", col_prefix = "Deblur.noPos_")

#OpenVnoPos
OpenVnoPos <- as.matrix(DM[grep("open.ref_", rownames(DM)), grep("Deblur.noPos", colnames(DM))])
OpenVnoPos_order <- reorder_pipeline_samples(OpenVnoPos, row_prefix = "open.ref_", col_prefix = "Deblur.noPos_")

###boxplots for DadaVDeblur, DadaVnoPos, UnoiseVDeblur, UnoiseVnoPos, OpenVDeblur, OpenVnoPos, DeblurVnoPos

noPOs_boxplot <- boxplot(diag(DavDe_ordered),
                         diag(DadaVnoPos_order),
                         diag(UvDe_ordered),
                         diag(UnoiseVnoPos_order),
                         diag(DeblurVOpen_order),
                         diag(OpenVnoPos_order),
                         names=c("DADA2_Deblur", "DADA2_noPos", "UNOISE3_Deblur", "UNOISE3_noPos", "Open_Deblur", "Open_noPos"),
                         ylab="Weighted UniFrac Distance")
noPos_weighted_unifrac <- recordPlot(noPOs_boxplot)

#Used to make figure 4
Unifrac_box <- recordPlot(Unifrac_boxplot)


###################################################################unweighted unifrac##############################################################
U_DM <- read.table("unweighted_unifrac_rare_Combined_Blueberry_fixed.txt", sep="\t", header=TRUE, row.names=1)
rownames(U_DM) <- gsub("-", ".", rownames(U_DM))

# Unoise vs Dada
U_unoiseVDada <- U_DM[grepl("Unoise*", names(U_DM)) ,]
U_unoiseVDada <- U_unoiseVDada[, grepl("Dada*", names(U_unoiseVDada))]
U_UvDa <- data.matrix(U_unoiseVDada)
# Already in order so don't need to re-order!

# Unoise vs Deblur
U_unoiseVDeblur <- U_DM[grepl("Unoise*", names(U_DM)),]
U_unoiseVDeblur <- U_unoiseVDeblur[, grepl("Deblur*", names(U_unoiseVDeblur))]
U_UvDe <- data.matrix(U_unoiseVDeblur)

U_UvDe_ordered <- reorder_pipeline_samples(dm_input=U_UvDe, row_prefix="Unoise_", col_prefix="Deblur_")

# Unoise vs Deblur
U_dadaVdeblur <- U_DM[grepl("Dada*", names(U_DM)) ,]
U_dadaVdeblur <- U_dadaVdeblur[,grepl("Deblur*", names(U_dadaVdeblur))]
U_DavDe <- data.matrix(U_dadaVdeblur)

U_DavDe_ordered <- reorder_pipeline_samples(dm_input=U_DavDe, row_prefix="Dada_", col_prefix="Deblur_")

# Dada vs Dada
U_dadaVdada <- U_DM[grep("Dada", names(U_DM)) ,grep("Dada", names(U_DM))]
U_DavDa <- data.matrix(U_dadaVdada)
U_DavDa_ordered <- reorder_pipeline_samples(dm_input=U_DavDa, row_prefix="Dada_", col_prefix="Dada_")

# Deblur vs Deblur
U_deblurVdeblur <- U_DM[grep("Deblur_", names(U_DM)) ,grep("Deblur_", names(U_DM))]
U_DevDe <- data.matrix(U_deblurVdeblur)
U_DevDe_ordered <- reorder_pipeline_samples(dm_input=U_DevDe, row_prefix="Deblur_", col_prefix="Deblur_")

# unoise vs unoise
U_unoiseVunoise <- U_DM[grep("Unoise", names(U_DM)) ,grep("Unoise", names(U_DM))]
U_UvU <- data.matrix(U_unoiseVunoise)
U_UvU_ordered <- reorder_pipeline_samples(dm_input=U_UvU, row_prefix="Unoise_", col_prefix="Unoise_")


#add in open
U_DadaVOpen <- U_DM[grep("open.ref_", names(U_DM)), grep("Dada_", names(U_DM))]
U_DadaVOpenorder <- reorder_pipeline_samples(dm_input=as.matrix(U_DadaVOpen), row_prefix="open.ref_", col_prefix ="Dada_" )

U_DeblurVOpen <- U_DM[grep("open.ref_", names(U_DM)), grep("Deblur_", names(U_DM))]
U_DeblurVOpen_order <- reorder_pipeline_samples(dm_input = as.matrix(U_DeblurVOpen), row_prefix = "open.ref_", col_prefix = "Deblur_")

U_UnoiseVOpen <- U_DM[grep("open.ref_", names(U_DM)), grep("Unoise_", names(U_DM))]
U_UnoiseVOpen_order <- reorder_pipeline_samples(dm_input = as.matrix(U_UnoiseVOpen), row_prefix = "open.ref_", col_prefix = "Unoise_")



# Sanity checks:
dim(U_UvDa)
dim(U_UvDe_ordered)
dim(U_DavDe_ordered)
dim(U_DadaVOpenorder)
dim(U_DeblurVOpen_order)
dim(U_UnoiseVOpen_order)

U_DM["Dada_Bact194", "Deblur_Bact194"]
U_DavDe_ordered["Dada_Bact194", "Deblur_Bact194"]
U_DM["Unoise_Bact1z2z2","Dada_Bact1z2z2"]
U_UvDa["Unoise_Bact1z2z2",which(rownames(U_UvDa) == "Unoise_Bact1z2z2")]

U_DM["Unoise_Bact1z2z2","Deblur_Bact1z2z2"]
U_UvDe_ordered["Unoise_Bact1z2z2",which(rownames(U_UvDe_ordered) == "Unoise_Bact1z2z2")]

#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))


U_Unifrac_boxplot <- boxplot(diag(U_DavDe_ordered),
                           diag(U_UvDa),
                           diag(U_UvDe_ordered),
                           diag(U_DadaVOpenorder),
                           diag(U_DeblurVOpen_order),
                           diag(U_UnoiseVOpen_order),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "DADA2_Open",  "Deblur_Open", "UNOISE3_Open"),
                           ylab="Unweighted UniFrac Distance",
                           cex.axis=0.66,
                           outline=FALSE)

#Used to make figure 4
blueberry_U_Unifrac_box <- recordPlot(U_Unifrac_boxplot)


#DeblurVDeblurnoPos
rownames(U_DM) <- gsub("001101", "", rownames(U_DM))
colnames(U_DM) <- gsub("001101", "", colnames(U_DM))
U_DeblurVDeblur.noPos <- as.matrix(U_DM[grep("Deblur_", rownames(U_DM)), grep("Deblur.noPos_",colnames(U_DM))])
U_DeblurVDeblur.noPos_order <- reorder_pipeline_samples(U_DeblurVDeblur.noPos, row_prefix = "Deblur_", col_prefix = "Deblur.noPos_")

#DadaVnoPOs 
U_DadaVnoPos <- as.matrix(U_DM[grep("Dada_", rownames(U_DM)), grep("Deblur.noPos_", colnames(U_DM))])
U_DadaVnoPos_order <- reorder_pipeline_samples(U_DadaVnoPos, row_prefix = "Dada_", col_prefix = "Deblur.noPos_")

#UnoiseVnoPos
U_UnoiseVnoPos <- as.matrix(U_DM[grep("Unoise_", rownames(U_DM)), grep("Deblur.noPos_", colnames(U_DM))])
U_UnoiseVnoPos_order <- reorder_pipeline_samples(U_UnoiseVnoPos, row_prefix = "Unoise_", col_prefix = "Deblur.noPos_")

#OpenVnoPos
U_OpenVnoPos <- as.matrix(U_DM[grep("open.ref_", rownames(U_DM)), grep("Deblur.noPos", colnames(U_DM))])
U_OpenVnoPos_order <- reorder_pipeline_samples(U_OpenVnoPos, row_prefix = "open.ref_", col_prefix = "Deblur.noPos_")

###boxplots for DadaVDeblur, DadaVnoPos, UnoiseVDeblur, UnoiseVnoPos, OpenVDeblur, OpenVnoPos, DeblurVnoPos

U_noPOs_boxplot <- boxplot(diag(U_DavDe_ordered),
                         diag(U_DadaVnoPos_order),
                         diag(U_UvDe_ordered),
                         diag(U_UnoiseVnoPos_order),
                         diag(U_DeblurVOpen_order),
                         diag(U_OpenVnoPos_order),
                         names=c("DADA2_Deblur", "DADA2_noPos", "UNOISE3_Deblur", "UNOISE3_noPos", "Open_Deblur", "Open_noPos"),
                         ylab="Unweighted UniFrac Distance")

noPos_unweighted <- recordPlot(U_noPOs_boxplot)
