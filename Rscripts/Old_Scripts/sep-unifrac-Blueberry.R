#Generate boxplots for weighted unifrac based on seperate DM tables

setwd("~/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/final_combined/sequences/pynast_aligned_seqs/sep_aligned/")

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

#read in all the distance matrixs
DadaVUnoise <- read.table("DadavUnoise/distances/weighted_unifrac_rare_DadavUnoise.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = F)
DadaVOpenref <- read.table("DadaVOpen-ref/distances/weighted_unifrac_rare_dadaVopen-ref.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
DeblurVDada <- read.table("DeblurVDada/distances/weighted_unifrac_rare_DeblurVDada.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
DeblurvNoPos <- read.table("DeblurVnoPos/distances/weighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
rownames(DeblurvNoPos) <- gsub("-", "", rownames(DeblurvNoPos))
colnames(DeblurvNoPos) <- gsub("-", "", rownames(DeblurvNoPos))
DeblurVopenref <- read.table("DeblurVOpen-ref/distances/weighted_unifrac_rare_deblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DeblurVUnoise <- read.table("DeblurVUnoise/distances/weighted_unifrac_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
UnoiseVOpenref <- read.table("UnoiseVOpen-ref/distances/weighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)


format_DM <- function(DM, row_prefix, col_prefix) {
  
  ret_tab <- DM[-grep(col_prefix, rownames(DM)), -grep(row_prefix, colnames(DM)),]
  ret_tab <- reorder_pipeline_samples(dm_input = ret_tab, row_prefix = row_prefix, col_prefix = col_prefix)
  
  return(ret_tab)
}

#prepare all the tables to be graphed.

DadaVUnoise_format <- as.matrix(format_DM(DM = DadaVUnoise, row_prefix = "Dada_", col_prefix = "Unoise_"))
dim(DadaVUnoise_format)
DadaVOpenref_format <- as.matrix(format_DM(DM=DadaVOpenref, row_prefix = "Dada_", col_prefix = "open.ref_"))
dim(DadaVOpenref_format)
DeblurVDada_format <- as.matrix(format_DM(DM=DeblurVDada, row_prefix="Deblur_", col_prefix = "Dada_"))
dim(DeblurVDada_format)
#fixing row and col names.
rownames(DeblurvNoPos) <- gsub("001101", "", rownames(DeblurvNoPos))
colnames(DeblurvNoPos) <- gsub("001101", "", colnames(DeblurvNoPos))
DeblurvNoPos_format <- as.matrix(format_DM(DM=DeblurvNoPos, row_prefix = "Deblur_", col_prefix = "DeblurnoPos_"))
DeblurVopenref_format <- as.matrix(format_DM(DM=DeblurVopenref, row_prefix = "Deblur_", col_prefix = "open-ref_"))
DeblurVUnoise_format <- as.matrix(format_DM(DM=DeblurVUnoise, row_prefix = "Deblur_", col_prefix = "Unoise_"))
UnoiseVOpenref_format <- as.matrix(format_DM(DM=UnoiseVOpenref, row_prefix = "Unoise_", col_prefix = "open-ref_"))





par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(DadaVUnoise_format),
                           diag(DeblurVDada_format),
                           diag(DeblurVUnoise_format),
                           diag(DeblurVopenref_format),
                           diag(UnoiseVOpenref_format),
                           diag(DadaVOpenref_format),
                           names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3","Deblur_Open-Ref", "Unoise_Open-Ref", "Dada_Open-Ref"),
                           ylab="Weighted UniFrac Distance",
                           cex.axis=0.5,
                           outline=FALSE)

Blueberry_weighted_boxplot <- recordPlot(Unifrac_boxplot)



####################################################### After low confidence filters

#read in all the distance matrixs
low_DadaVUnoise <- read.table("DadavUnoise/low_conf_filt/distances/weighted_unifrac_filt_rare_DadavUnoise.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = F)
low_DadaVOpenref <- read.table("DadaVOpen-ref/low_conf_filt/distances/weighted_unifrac_filt_rare_dadaVopen-ref.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
low_DeblurVDada <- read.table("DeblurVDada/low_conf_filt/distances/weighted_unifrac_filt_rare_DeblurVDada.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
low_DeblurvNoPos <- read.table("DeblurVnoPos/low_conf_filt/distances/weighted_unifrac_filt_rare_DEblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
rownames(low_DeblurvNoPos) <- gsub("-", "", rownames(low_DeblurvNoPos))
colnames(low_DeblurvNoPos) <- gsub("-", "", rownames(low_DeblurvNoPos))
low_DeblurVopenref <- read.table("DeblurVOpen-ref/low_conf_filt/distances/weighted_unifrac_filt_rare_deblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
low_DeblurVUnoise <- read.table("DeblurVUnoise/low_conf_filt/distances/weighted_unifrac_filt_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
low_UnoiseVOpenref <- read.table("UnoiseVOpen-ref/low_conf_filt/distances/weighted_unifrac_filt_rare_UnoiseVOpen_ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)



###format the tables
low_DadaVUnoise_format <- as.matrix(format_DM(DM = low_DadaVUnoise, row_prefix = "Dada_", col_prefix = "Unoise_"))
dim(low_DadaVUnoise_format)
low_DadaVOpenref_format <- as.matrix(format_DM(DM=low_DadaVOpenref, row_prefix = "Dada_", col_prefix = "open.ref_"))
dim(low_DadaVOpenref_format)
low_DeblurVDada_format <- as.matrix(format_DM(DM=low_DeblurVDada, row_prefix="Deblur_", col_prefix = "Dada_"))
dim(low_DeblurVDada_format)
#fixing row and col names.
rownames(low_DeblurvNoPos) <- gsub("001101", "", rownames(low_DeblurvNoPos))
colnames(low_DeblurvNoPos) <- gsub("001101", "", colnames(low_DeblurvNoPos))
low_DeblurvNoPos_format <- as.matrix(format_DM(DM=low_DeblurvNoPos, row_prefix = "Deblur_", col_prefix = "DeblurnoPos_"))
low_DeblurVopenref_format <- as.matrix(format_DM(DM=low_DeblurVopenref, row_prefix = "Deblur_", col_prefix = "open-ref_"))
low_DeblurVUnoise_format <- as.matrix(format_DM(DM=low_DeblurVUnoise, row_prefix = "Deblur_", col_prefix = "Unoise_"))
low_UnoiseVOpenref_format <- as.matrix(format_DM(DM=low_UnoiseVOpenref, row_prefix = "Unoise_", col_prefix = "open-ref_"))



###################################################graph it

Unifrac_boxplot <- boxplot(diag(low_DadaVUnoise_format),
                           diag(low_DeblurVDada_format),
                           diag(low_DeblurVUnoise_format),
                           diag(low_DeblurvNoPos_format),
                           diag(low_DeblurVopenref_format),
                           diag(low_UnoiseVOpenref_format),
                           diag(low_DadaVOpenref_format),
                           names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_Deb_noPos", "Deblur_Open", "Unoise_Open", "Dada_Open"),
                           ylab="Weighted UniFrac Distance", 
                           outline=FALSE)


########################################### Unweighted Unifrac ##################################################################

#read in all the distance matrixs
U_DadaVUnoise <- read.table("DadavUnoise/distances/unweighted_unifrac_rare_DadavUnoise.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = F)
U_DadaVOpenref <- read.table("DadaVOpen-ref/distances/unweighted_unifrac_rare_dadaVopen-ref.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
U_DeblurVDada <- read.table("DeblurVDada/distances/unweighted_unifrac_rare_DeblurVDada.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F)
U_DeblurvNoPos <- read.table("DeblurVnoPos/distances/unweighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
rownames(U_DeblurvNoPos) <- gsub("-", "", rownames(U_DeblurvNoPos))
colnames(U_DeblurvNoPos) <- gsub("-", "", rownames(U_DeblurvNoPos))
U_DeblurVopenref <- read.table("DeblurVOpen-ref/distances/unweighted_unifrac_rare_deblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DeblurVUnoise <- read.table("DeblurVUnoise/distances/unweighted_unifrac_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_UnoiseVOpenref <- read.table("UnoiseVOpen-ref/distances/unweighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)



#format the tables

U_DadaVUnoise_format <- as.matrix(format_DM(DM = U_DadaVUnoise, row_prefix = "Dada_", col_prefix = "Unoise_"))
dim(U_DadaVUnoise_format)
U_DadaVOpenref_format <- as.matrix(format_DM(DM=U_DadaVOpenref, row_prefix = "Dada_", col_prefix = "open.ref_"))
dim(U_DadaVOpenref_format)
U_DeblurVDada_format <- as.matrix(format_DM(DM=U_DeblurVDada, row_prefix="Deblur_", col_prefix = "Dada_"))
dim(U_DeblurVDada_format)
#fixing row and col names.
rownames(U_DeblurvNoPos) <- gsub("001101", "", rownames(U_DeblurvNoPos))
colnames(U_DeblurvNoPos) <- gsub("001101", "", colnames(U_DeblurvNoPos))
U_DeblurvNoPos_format <- as.matrix(format_DM(DM=U_DeblurvNoPos, row_prefix = "Deblur_", col_prefix = "DeblurnoPos_"))
U_DeblurVopenref_format <- as.matrix(format_DM(DM=U_DeblurVopenref, row_prefix = "Deblur_", col_prefix = "open-ref_"))
U_DeblurVUnoise_format <- as.matrix(format_DM(DM=U_DeblurVUnoise, row_prefix = "Deblur_", col_prefix = "Unoise_"))
U_UnoiseVOpenref_format <- as.matrix(format_DM(DM=U_UnoiseVOpenref, row_prefix = "Unoise_", col_prefix = "open-ref_"))



#### graph it

Unifrac_boxplot <- boxplot(diag(U_DadaVUnoise_format),
                           diag(U_DeblurVDada_format),
                           diag(U_DeblurVUnoise_format),
                           diag(U_DeblurvNoPos_format),
                           diag(U_DeblurVopenref_format),
                           diag(U_UnoiseVOpenref_format),
                           diag(U_DadaVOpenref_format),
                           names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_Deb_noPos", "Deblur_Open", "Unoise_Open", "Dada_Open"),
                           ylab="Unweighted UniFrac Distance", 
                           outline=FALSE)


#format distance matrix for mantel tests on weighted unifrac

#testing sample to sample differences
