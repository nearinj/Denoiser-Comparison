#unweighted_unifrac_biscuit distances
####################################################### new dm's that include open-ref and no postive 


setwd("~/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/final_combined/distances/")


# Read in matrix of weighted unifrac distances at ASV level for samples across all 3 pipelines.
new_DM <- read.table("unweighted_unifrac_rare_CombinedV2_biscuit_taxa.txt", sep="\t", header=TRUE, row.names=1)

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
new_unoiseVDada <- new_DM[grepl("Unoise*", names(new_DM)) ,]
new_unoiseVDada <- new_unoiseVDada[, grepl("Dada*", names(new_unoiseVDada))]
new_UvDa <- data.matrix(new_unoiseVDada)
# Already in order so don't need to re-order!
dim(new_UvDa)
# Unoise vs Deblur
new_unoiseVDeblur <- new_DM[grepl("Unoise*", names(new_DM)),]
new_unoiseVDeblur <- new_unoiseVDeblur[, grepl("Deblur_E", names(new_unoiseVDeblur))]
new_UvDe <- data.matrix(new_unoiseVDeblur)
dim(new_UvDe)

new_UvDe_ordered <- reorder_pipeline_samples(dm_input=new_UvDe, row_prefix="Unoise__", col_prefix="Deblur_")

# Unoise vs Deblur
new_dadaVdeblur <- new_DM[grepl("Dada*", names(new_DM)) ,]
new_dadaVdeblur <- new_dadaVdeblur[,grepl("Deblur_E", names(new_dadaVdeblur))]
new_DavDe <- data.matrix(new_dadaVdeblur)

new_DavDe_ordered <- reorder_pipeline_samples(dm_input=new_DavDe, row_prefix="Dada_", col_prefix="Deblur_")

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
dim(new_UvDa)
dim(new_UvDe_ordered)
dim(new_DavDe_ordered)

gsub("Dada_", "", colnames(new_UvDa)) == gsub("Unoise__", "", rownames(new_UvDa))
gsub("Deblur_", "", colnames(new_UvDe_ordered)) == gsub("Unoise__", "", rownames(new_UvDe_ordered))
gsub("Deblur_", "", colnames(new_DavDe_ordered)) == gsub("Dada_", "", rownames(new_DavDe_ordered))



colnames(new_DM) <- gsub("001101", "", colnames(new_DM))
colnames(new_DM)
#postive deblur vs no postive
new_noPDeblurVdeblur <- new_DM[grepl("Deblur.no", names(new_DM)) ,]
rownames(new_noPDeblurVdeblur) <- gsub("001101", "", rownames(new_noPDeblurVdeblur))
rownames(new_noPDeblurVdeblur)
new_noPDeblurVdeblur <- new_noPDeblurVdeblur[,grepl("Deblur_", names(new_noPDeblurVdeblur))]
colnames(new_noPDeblurVdeblur)
new_noPDevDe <- data.matrix(new_noPDeblurVdeblur)
dim(new_noPDevDe)

new_noPDevDe_ordered <- reorder_pipeline_samples(dm_input=new_noPDevDe, row_prefix="Deblur-noPos_", col_prefix="Deblur_")

#no postive vs DADA
new_noPDeblurVDADA <- new_DM[grepl("Deblur.no", names(new_DM)) ,]
rownames(new_noPDeblurVDADA) <- gsub("001101", "", rownames(new_noPDeblurVDADA))
rownames(new_noPDeblurVDADA)
new_noPDeblurVDADA <- new_noPDeblurVDADA[,grepl("Dada_", names(new_noPDeblurVDADA))]
colnames(new_noPDeblurVDADA)
new_noPDevDa <- data.matrix(new_noPDeblurVDADA)
dim(new_noPDevDa)

new_noPDevDa_ordered <- reorder_pipeline_samples(dm_input=new_noPDevDa, row_prefix="Deblur-noPos_", col_prefix="Dada_")


#no postive vs UNOISE
new_noPDeblurVUNOISE <- new_DM[grepl("Deblur.no", names(new_DM)) ,]
rownames(new_noPDeblurVUNOISE) <- gsub("001101", "", rownames(new_noPDeblurVUNOISE))
rownames(new_noPDeblurVUNOISE)
new_noPDeblurVUNOISE <- new_noPDeblurVUNOISE[,grepl("Unoise_", names(new_noPDeblurVUNOISE))]
colnames(new_noPDeblurVUNOISE)
new_noPDevU <- data.matrix(new_noPDeblurVUNOISE)
dim(new_noPDevU)

new_noPDevU_ordered <- reorder_pipeline_samples(dm_input=new_noPDevU, row_prefix="Deblur-noPos_", col_prefix="Unoise__")


#open-ref v DADA2

new_openRefvDADA <- new_DM[grepl("open.ref_", names(new_DM)) ,]
rownames(new_openRefvDADA)
new_openRefvDADA <- new_openRefvDADA[,grepl("Dada_", names(new_openRefvDADA))]
colnames(new_openRefvDADA)
#remove extra samples that passed the filtering ... this was previously done manually.
new_openRefvDADA <- new_openRefvDADA[gsub("Dada_", "open-ref_", colnames(new_openRefvDADA)),]
rownames(new_openRefvDADA)
new_OPRvDA <- data.matrix(new_openRefvDADA)
dim(new_OPRvDA)

new_OPRvDA_ordered <- reorder_pipeline_samples(dm_input=new_OPRvDA, row_prefix = "open-ref_", col_prefix = "Dada_")
colnames(new_OPRvDA_ordered)

#sanity check
which(gsub("open-ref_", "", rownames(new_OPRvDA_ordered)) == gsub("Dada_", "", colnames(new_OPRvDA_ordered)))


#open-ref v Deblur
new_openRefvDeblur <- new_DM[grepl("open.ref_", names(new_DM)),]
rownames(new_openRefvDeblur)
new_openRefvDeblur<- new_openRefvDeblur[,grepl("Deblur_", names(new_openRefvDeblur))]
colnames(new_openRefvDeblur)
#remove extra samples that passed the filtering rarefaction..
new_openRefvDeblur <- new_openRefvDeblur[gsub("Deblur_", "open-ref_", colnames(new_openRefvDeblur)),]
rownames(new_openRefvDeblur)
new_OPRvDE <- data.matrix(new_openRefvDeblur)
dim(new_OPRvDE)


new_OPRvDE_ordered <- reorder_pipeline_samples(dm_input=new_OPRvDE, row_prefix = "open-ref_", col_prefix = "Deblur_")
colnames(new_OPRvDE_ordered)

#sanity check
which(gsub("open-ref_", "", rownames(new_OPRvDE_ordered)) == gsub("Deblur_", "", colnames(new_OPRvDE_ordered)))


#open-ref v Unoise
new_openRefvUNOISE <- new_DM[grepl("open.ref_", names(new_DM)),]
rownames(new_openRefvUNOISE)
new_openRefvUNOISE <- new_openRefvUNOISE[,grepl("Unoise_", names(new_openRefvUNOISE))]
colnames(new_openRefvUNOISE)
#remove extra samples that passed the filtering rarefaction
new_openRefvUNOISE <- new_openRefvUNOISE[gsub("Unoise__", "open-ref_", colnames(new_openRefvUNOISE)),]
rownames(new_openRefvUNOISE)
new_OPRvUN <- data.matrix(new_openRefvUNOISE)
dim(new_OPRvUN)

new_OPRvUN_ordered <- reorder_pipeline_samples(dm_input = new_OPRvUN, row_prefix = "open-ref_", col_prefix="Unoise__")

which(gsub("open-ref_", "", rownames(new_OPRvUN_ordered)) == gsub("Unoise__", "", colnames(new_OPRvUN_ordered)))




#sanity checks
new_DM["Unoise__ERR2042060", "open.ref_ERR2042060"] == new_OPRvUN_ordered["open-ref_ERR2042060", "Unoise__ERR2042060"]

new_OPRvUN_ordered["open-ref_ERR2042060", "Unoise__ERR2042060"]
rownames(new_OPRvUN_ordered)
#box plot of intre samples distances
par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(new_DavDe_ordered),
                           diag(new_UvDa),
                           diag(new_UvDe_ordered), 
                           diag(new_noPDevDe_ordered),
                           diag(new_noPDevDa_ordered),
                           diag(new_noPDevU_ordered),
                           diag(new_OPRvDA),
                           diag(new_OPRvDE),
                           diag(new_OPRvUN),
                           names=c("DADA2_Deblur", "DADA2_UNOISE3", "UNOISE3_Deblur", "NoPosDeb_Deblur", "NoPosDeb_DADA2", "NoPosDeb_UNOISE3", "Open_DADA2", 
                                   "Open_Deblur", "Open_UNOISE3"),
                           ylab="Unweighted UniFrac Distance", 
                           outline=FALSE)
# distances are not exactly the same.... why is this???? (for the orginal three...)
# distance difference between postive and non postive not that large.... 
# look into this tomorrow.