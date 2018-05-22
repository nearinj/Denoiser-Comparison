#sep-unifrac for BISCUIT data

#functions

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

format_DM <- function(DM, row_prefix, col_prefix) {
  
  ret_tab <- DM[-grep(col_prefix, rownames(DM)), -grep(row_prefix, colnames(DM)),]
  ret_tab <- reorder_pipeline_samples(dm_input = ret_tab, row_prefix = row_prefix, col_prefix = col_prefix)
  
  return(ret_tab)
}

#read in all the tables

setwd("~/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/final_combined/sequences/pynast_aligned_seqs/sep-align/")

DadaVDeblur <- read.table("DadaVDeblur/distances/weighted_unifrac_rare_DadaVDeblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DadaVUnoise <- read.table("DadaVUnoise/distances/weighted_unifrac_rare_DadaVUnoise.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
DadaVOpenRef <- read.table("DadaVOpen-ref/distances/weighted_unifrac_rare_DadavOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
UnoiseVOpenRef <- read.table("UnoiseVOpen-ref/distances/weighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DeblurVOpenRef <- read.table("DeblurVOpen-ref/distances/weighted_unifrac_rare_DeblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DeblurVnoPos <- read.table("DeblurVnoPos/distances/weighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
colnames(DeblurVnoPos) <- gsub("001101", "", colnames(DeblurVnoPos))
rownames(DeblurVnoPos) <- gsub("001101", "", rownames(DeblurVnoPos))
DeblurVUnoise <- read.table("DeblurVUnoise/distances/weighted_unifrac_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)

DadaVDeblur_format <- as.matrix(format_DM(DadaVDeblur, "Deblur_", "Dada_"))
DadaVUnoise_format <- as.matrix(format_DM(DadaVUnoise, "Dada_", "Unoise__"))
DadaVOpenRef_format <- as.matrix(format_DM(DadaVOpenRef, "open-ref_", "Dada_"))
UnoiseVOpenRef_format <- as.matrix(format_DM(UnoiseVOpenRef, "open-ref_", "Unoise__"))
DeblurVOpenRef_format <- as.matrix(format_DM(DeblurVOpenRef, "Deblur_", "open-ref_"))
DeblurVnoPos_format <- as.matrix(format_DM(DeblurVnoPos, "Deblur-noPos_", "Deblur_"))
DeblurVUnoise_format <- as.matrix(format_DM(DeblurVUnoise, "Deblur_", "Unoise__"))


weighted_biscuit <- boxplot(diag(DadaVUnoise_format), 
        diag(DadaVDeblur_format),
        diag(DeblurVUnoise_format),
        diag(DeblurVOpenRef_format),
        diag(UnoiseVOpenRef_format),
        diag(DadaVOpenRef_format),
        diag(DeblurVnoPos_format),
        names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_Open", "Unoise_Open", "Dada_Open", "Deblur_noPos"),
        ylab="Weighted UniFrac Distance",
        cex.axis=0.5,
        outline=FALSE)

recorded_weighted_biscuit <- recordPlot(weighted_biscuit)

######################## Unweighted Unifrac

#read in unweighted tables

U_DadaVDeblur <- read.table("DadaVDeblur/distances/unweighted_unifrac_rare_DadaVDeblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DadaVUnoise <- read.table("DadaVUnoise/distances/unweighted_unifrac_rare_DadaVUnoise.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
U_DadaVOpenRef <- read.table("DadaVOpen-ref/distances/unweighted_unifrac_rare_DadavOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_UnoiseVOpenRef <- read.table("UnoiseVOpen-ref/distances/unweighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DeblurVOpenRef <- read.table("DeblurVOpen-ref/distances/unweighted_unifrac_rare_DeblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DeblurVnoPos <- read.table("DeblurVnoPos/distances/unweighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
colnames(U_DeblurVnoPos) <- gsub("001101", "", colnames(U_DeblurVnoPos))
rownames(U_DeblurVnoPos) <- gsub("001101", "", rownames(U_DeblurVnoPos))
U_DeblurVUnoise <- read.table("DeblurVUnoise/distances/unweighted_unifrac_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)

#format the unweighted tables

U_DadaVDeblur_format <- as.matrix(format_DM(U_DadaVDeblur, "Deblur_", "Dada_"))
U_DadaVUnoise_format <- as.matrix(format_DM(U_DadaVUnoise, "Dada_", "Unoise__"))
U_DadaVOpenRef_format <- as.matrix(format_DM(U_DadaVOpenRef, "open-ref_", "Dada_"))
U_UnoiseVOpenRef_format <- as.matrix(format_DM(U_UnoiseVOpenRef, "open-ref_", "Unoise__"))
U_DeblurVOpenRef_format <- as.matrix(format_DM(U_DeblurVOpenRef, "Deblur_", "open-ref_"))
U_DeblurVnoPos_format <- as.matrix(format_DM(U_DeblurVnoPos, "Deblur-noPos_", "Deblur_"))
U_DeblurVUnoise_format <- as.matrix(format_DM(U_DeblurVUnoise, "Deblur_", "Unoise__"))

#graph it

boxplot(diag(U_DadaVUnoise_format), 
        diag(U_DadaVDeblur_format),
        diag(U_DeblurVUnoise_format),
        diag(U_DeblurVnoPos_format),
        diag(U_DeblurVOpenRef_format),
        diag(U_UnoiseVOpenRef_format),
        diag(U_DadaVOpenRef_format),
        names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_Deb_noPos", "Deblur_Open", "Unoise_Open", "Dada_Open"),
        ylab="Unweighted UniFrac Distance", 
        outline=FALSE)
