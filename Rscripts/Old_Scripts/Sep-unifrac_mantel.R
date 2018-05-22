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



#shuffle sequence trouble shooting
DadaVOpenref_shuffle <- read.table("DadaVOpen-ref/test/distances/weighted_unifrac_DadaVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DadaVOpenref_shuffle_format <- as.matrix(format_DM(DM=DadaVOpenref_shuffle, row_prefix = "Dada_", col_prefix = "open-ref_"))
shuffle_plot <- boxplot(diag(DadaVOpenref_format), diag(DadaVOpenref_shuffle_format))

par(xpd= NA, bg = "transparent", oma= c(2,2,0,0))
Unifrac_boxplot <- boxplot(diag(DadaVUnoise_format),
                           diag(DeblurVDada_format),
                           diag(DeblurVUnoise_format),
                           diag(DeblurvNoPos_format),
                           diag(DeblurVopenref_format),
                           diag(UnoiseVOpenref_format),
                           diag(DadaVOpenref_format),
                           names=c("Dada2_UNOISE3", "Deblur_Dada2", "Deblur_UNOISE3", "Deblur_Deb_noPos", "Deblur_Open", "Unoise_Open", "Dada_Open"),
                           ylab="Weighted UniFrac Distance", 
                           outline=FALSE)




diag(DadavUnoise_format)
DadavUnoise_format



#format distance matrix for mantel tests
library(vegan)

mantel_format <- function(tab, prefix){
  
  ret_tab <- tab[grepl(prefix, rownames(tab)), grepl(prefix, colnames(tab))]
  
  return(ret_tab)
}

mantel_test <- function(tab, prefix1, prefix2){
  prefix1_tab <- mantel_format(tab, prefix1)
  prefix2_tab <- mantel_format(tab, prefix2)
  prefix1_tab <- prefix1_tab[gsub(prefix2, prefix1, rownames(prefix2_tab)), gsub(prefix2, prefix1, colnames(prefix2_tab))]
  #check if they are the same
  if(length(which(gsub(prefix1, "", rownames(prefix1_tab)) != gsub(prefix2,"", rownames(prefix2_tab)))) == length(rownames(prefix1_tab))){
    return("Error")
  }else if(length(which(gsub(prefix1,"", colnames(prefix1_tab)) != gsub(prefix2,"", colnames(prefix2_tab)))) == length(colnames(prefix2_tab))){
    return("Error")
  }else
    return(mantel(prefix1_tab, prefix2_tab, permutations = 999, method="pearson"))
  
}
Mantel_DadaVOpenref_Dada <- mantel_format(DadaVOpenref, "Dada") 
Mantel_DadaVOpenref_Open <- mantel_format(DadaVOpenref, "open")
#check if they are in the same order
Mantel_DadaVOpenref_Open <- Mantel_DadaVOpenref_Open[gsub("Dada_", "open-ref_", rownames(Mantel_DadaVOpenref_Dada)), gsub("Dada_", "open.ref_", colnames(Mantel_DadaVOpenref_Dada))]
gsub("Dada_","",rownames(Mantel_DadaVOpenref_Dada)) == gsub("open-ref_", "", rownames(Mantel_DadaVOpenref_Open))
gsub("Dada_", "", colnames(Mantel_DadaVOpenref_Dada)) == gsub("open.ref_", "", colnames(Mantel_DadaVOpenref_Open))


DadaVopenRef_Mantel_stats <- mantel(Mantel_DadaVOpenref_Dada, Mantel_DadaVOpenref_Open, method="pearson", permutations=999)
DadaVopenRef_Mantel_stats

DadaVUnoise_mantel <- mantel_test(DadaVUnoise, "Dada_", "Unoise_")
DadaVUnoise_mantel

DeblurVDada_mantel <- mantel_test(DeblurVDada, "Deblur_", "Dada_")
DeblurVDada_mantel

DeblurvNoPos_mantel <- mantel_test(DeblurvNoPos, "Deblur_", "DeblurnoPos_")
DeblurvNoPos_mantel

DeblurVopenref_mantel <- mantel_test(DeblurVopenref, "Deblur_", "open-ref_")
DeblurVopenref_mantel

DeblurVUnoise_mantel <- mantel_test(DeblurVUnoise, "Deblur_", "Unoise_") 
DeblurVUnoise_mantel

UnoiseVOpenref_mantel <- mantel_test(UnoiseVOpenref, "Unoise_", "open-ref_")
UnoiseVOpenref_mantel

mantel_r2_table <- data.frame(Dada=double(4),
                              Deblur=double(4),
                              Unoise=double(4),
                              Openref=double(4)
                              )
rownames(mantel_r2_table) <- c("Dada", "Unoise", "Deblur", "Openref")

mantel_r2_table["Dada", "Unoise"] <- DadaVUnoise_mantel$statistic*DadaVUnoise_mantel$statistic
mantel_r2_table["Dada", "Deblur"] <- DeblurVDada_mantel$statistic*DeblurVDada_mantel$statistic
mantel_r2_table["Dada", "Openref"] <- DadaVopenRef_Mantel_stats$statistic <- DadaVopenRef_Mantel_stats$statistic
mantel_r2_table["Unoise", "Dada"] <- DadaVUnoise_mantel$statistic*DadaVUnoise_mantel$statistic
mantel_r2_table["Unoise", "Deblur"] <- DeblurVUnoise_mantel$statistic*DeblurVUnoise_mantel$statistic
mantel_r2_table["Unoise", "Openref"] <- DeblurVopenref_mantel$statistic*DeblurVopenref_mantel$statistic
mantel_r2_table["Deblur", "Dada"] <- DeblurVDada_mantel$statistic*DeblurVDada_mantel$statistic
mantel_r2_table["Deblur", "Unoise"] <- DeblurVUnoise_mantel$statistic*DeblurVUnoise_mantel$statistic
mantel_r2_table["Deblur", "Openref"] <- DeblurVopenref_mantel$statistic*DeblurVopenref_mantel$statistic
mantel_r2_table["Openref", "Deblur"] <- DeblurVopenref_mantel$statistic*DeblurVopenref_mantel$statistic
mantel_r2_table["Openref", "Unoise"] <- DeblurVopenref_mantel$statistic*DeblurVopenref_mantel$statistic
mantel_r2_table["Openref", "Dada"] <-  DadaVopenRef_Mantel_stats$statistic <- DadaVopenRef_Mantel_stats$statistic
mantel_r2_table <- mantel_r2_table[rownames(mantel_r2_table), rownames(mantel_r2_table)]
mantel_r2_matrix <- as.matrix(mantel_r2_table)

library(gplots)
heatmap.2(mantel_r2_matrix, dendrogram = "none")


