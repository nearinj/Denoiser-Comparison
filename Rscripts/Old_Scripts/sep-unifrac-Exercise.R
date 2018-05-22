#sep-unifrac Exercise

#then im done with this bullshit.


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

#read in the tables

setwd("~/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/final_combined/sequences/pynast_aligned_seqs/sep-align/")



DadaVDeblur <- read.table("DadaVDeblur/distances/weighted_unifrac_rare_DadaVDeblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DadaVUnoise <- read.table("DadaVUnoise/distances/weighted_unifrac_rare_DadaVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DadaVOpenref <- read.table("DadaVOpen-ref/distances/weighted_unifrac_rare_DadaVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DeblurVUnoise <- read.table("UnoiseVDeblur/distances/weighted_unifrac_rare_UnoiseVDeblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DeblurVnoPos <- read.table("DeblurVnoPos/distances/weighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
colnames(DeblurVnoPos) <- gsub("001101", "", colnames(DeblurVnoPos))
rownames(DeblurVnoPos) <- gsub("001101", "", rownames(DeblurVnoPos))
DeblurVOpenref <- read.table("DeblurVOpen-ref/distances/weighted_unifrac_rare_DeblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
UnoiseVOpenref <- read.table("UnoiseVOpen-ref/distances/weighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)

#format the tables for the boxplot

DadaVDeblur_format <- as.matrix(format_DM(DadaVDeblur, "Deblur_", "Dada_"))
DadaVUnoise_format <- as.matrix(format_DM(DadaVUnoise, "Dada_", "Unoise_"))
DadaVOpenref_format <- as.matrix(format_DM(DadaVOpenref, "open-ref_", "Dada_"))
DeblurVUnoise_format <- as.matrix(format_DM(DeblurVUnoise, "Deblur_", "Unoise_"))
DeblurVnoPos_format <- as.matrix(format_DM(DeblurVnoPos, "Deblur-noPos_", "Deblur_"))
DeblurVOpenref_format <- as.matrix(format_DM(DeblurVOpenref, "Deblur_", "open-ref_"))
UnoiseVOpenref_format <- as.matrix(format_DM(UnoiseVOpenref, "open-ref_", "Unoise_"))



#graph it

Exercise_weighted <- boxplot(diag(DadaVUnoise_format), 
        diag(DadaVDeblur_format),
        diag(DeblurVUnoise_format),
        diag(DeblurVOpenref_format),
        diag(UnoiseVOpenref_format),
        diag(DadaVOpenref_format),
        diag(DeblurVnoPos_format),
        names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_Open", "Unoise_Open", "Dada_Open", "Deblur_noPos"),
        ylab="Weighted UniFrac Distance",
        cex.axis=0.5,
        outline=FALSE)


Exercise_weighted_boxplot <- recordPlot(Exercise_weighted)
###################################### unweighted unifrac

U_DadaVDeblur <- read.table("DadaVDeblur/distances/unweighted_unifrac_rare_DadaVDeblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DadaVUnoise <- read.table("DadaVUnoise/distances/unweighted_unifrac_rare_DadaVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DadaVOpenref <- read.table("DadaVOpen-ref/distances/unweighted_unifrac_rare_DadaVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DeblurVUnoise <- read.table("UnoiseVDeblur/distances/unweighted_unifrac_rare_UnoiseVDeblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DeblurVnoPos <- read.table("DeblurVnoPos/distances/unweighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
colnames(U_DeblurVnoPos) <- gsub("001101", "", colnames(U_DeblurVnoPos))
rownames(U_DeblurVnoPos) <- gsub("001101", "", rownames(U_DeblurVnoPos))
U_DeblurVOpenref <- read.table("DeblurVOpen-ref/distances/unweighted_unifrac_rare_DeblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_UnoiseVOpenref <- read.table("UnoiseVOpen-ref/distances/unweighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)

######################### format table

U_DadaVDeblur_format <- as.matrix(format_DM(U_DadaVDeblur, "Deblur_", "Dada_"))
U_DadaVUnoise_format <- as.matrix(format_DM(U_DadaVUnoise, "Dada_", "Unoise_"))
U_DadaVOpenref_format <- as.matrix(format_DM(U_DadaVOpenref, "open-ref_", "Dada_"))
U_DeblurVUnoise_format <- as.matrix(format_DM(U_DeblurVUnoise, "Deblur_", "Unoise_"))
U_DeblurVnoPos_format <- as.matrix(format_DM(U_DeblurVnoPos, "Deblur-noPos_", "Deblur_"))
U_DeblurVOpenref_format <- as.matrix(format_DM(U_DeblurVOpenref, "Deblur_", "open-ref_"))
U_UnoiseVOpenref_format <- as.matrix(format_DM(U_UnoiseVOpenref, "open-ref_", "Unoise_"))

#graph it

boxplot(diag(U_DadaVUnoise_format), 
        diag(U_DadaVDeblur_format),
        diag(U_DeblurVUnoise_format),
        diag(U_DeblurVnoPos_format),
        diag(U_DeblurVOpenref_format),
        diag(U_UnoiseVOpenref_format),
        diag(U_DadaVOpenref_format),
        names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_noPos", "Deblur_Open", "Unoise_Open", "Dada_Open"),
        ylab="Unweighted UniFrac Distance",
        outline=FALSE)
