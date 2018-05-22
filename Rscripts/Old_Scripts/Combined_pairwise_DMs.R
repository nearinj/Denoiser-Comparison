#make final distance matrix for PCoA plot

#functions
library(dplyr)
library(reshape2)
library(optparse)


format_DM <- function(DM, row_prefix, col_prefix) {
  
  ret_tab <- DM[-grep(col_prefix, rownames(DM)), -grep(row_prefix, colnames(DM)),]
  ret_tab <- reorder_pipeline_samples(dm_input = ret_tab, row_prefix = row_prefix, col_prefix = col_prefix)
  
  return(ret_tab)
}

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

#converts DM to pairwise table
#currently only works for same row and col names
to_pair_wise <- function(tab, not_same){
  #need to change this to get the combinations from the row and col names
  ret_tab <- data.frame(expand.grid(rownames(tab),colnames(tab)), stringsAsFactors = F)
  ret_tab$distance <- unlist(tab)
  colnames(ret_tab) <- c("S1", "S2", "distance")
  
  pass <- T
  #check if they all still match
  for(row in 1:nrow(ret_tab)){
    if(ret_tab[row,"distance"] == tab[as.character(ret_tab[row,"S1"]),as.character(ret_tab[row,"S2"])]){
    }else{
      break()	
      pass <- F
      print("failure")
    }
  }
  #make other pairings
  if(not_same){
    ret_tab2 <- ret_tab
    colnames(ret_tab2) <- c("S2", "S1", "distance")
    ret_tab <- rbind(ret_tab, ret_tab2)
  }
  
  if(pass){return(ret_tab)}
  else(return("fail"))
}

setwd("~/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/final_combined/sequences/pynast_aligned_seqs/sep_aligned/")

#load in weighted tables

DadaVUnoise <- read.table("DadavUnoise/distances/weighted_unifrac_rare_DadavUnoise.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = F, check.names=F)
DadaVOpenref <- read.table("DadaVOpen-ref/distances/weighted_unifrac_rare_dadaVopen-ref.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F, check.names=F)
DeblurVDada <- read.table("DeblurVDada/distances/weighted_unifrac_rare_DeblurVDada.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F, check.names=F)
DeblurvNoPos <- read.table("DeblurVnoPos/distances/weighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
colnames(DeblurvNoPos) <- gsub("001101", "", colnames(DeblurvNoPos))
rownames(DeblurvNoPos) <- gsub("001101", "", rownames(DeblurvNoPos))
DeblurVopenref <- read.table("DeblurVOpen-ref/distances/weighted_unifrac_rare_deblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
DeblurVUnoise <- read.table("DeblurVUnoise/distances/weighted_unifrac_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
UnoiseVOpenref <- read.table("UnoiseVOpen-ref/distances/weighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
Dada_only <- read.table("Dada/distances/weighted_unifrac_dada2_header.txt", sep="\t", header=T, row.names = 1, check.names = F, stringsAsFactors = F)
Unoise_only <- read.table("Unoise/distances/weighted_unifrac_unoise_header.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
Deblur_only <- read.table("Deblur/distances/weighted_unifrac_deblur_header.txt", sep="\t", header=T, row.names = 1, check.names = F, stringsAsFactors = F)
OpenRef_only <- read.table("Open/distances/weighted_unifrac_rare_Open-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
noPos_only <- read.table("noPos/distances/weighted_unifrac_rare_Deblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)

DadaVUnoise_format <- format_DM(DadaVUnoise, "Dada_", "Unoise_")
DadaVOpenref_format <- format_DM(DadaVOpenref, "Dada_", "open-ref_")
DeblurVDada_format <- format_DM(DeblurVDada, "Dada_", "Deblur_")
DeblurVNoPos_format <- format_DM(DeblurvNoPos, "Deblur_", "Deblur-noPos_")
DeblurVopenref_format <- format_DM(DeblurVopenref, "Deblur_", "open-ref_")
DeblurVUnoise_format <- format_DM(DeblurVUnoise, "Deblur_", "Unoise_")
UnoiseVOpenref_format <- format_DM(UnoiseVOpenref, "Unoise_", "open-ref_")

sample_name_list <- gsub("Deblur_", "", colnames(Deblur_only))


Dada_only_pairs <- to_pair_wise(Dada_only, F)
Unoise_only_pairs <- to_pair_wise(Unoise_only, F)
Deblur_only_pairs <- to_pair_wise(Deblur_only, F)
noPos_only_pairs <- to_pair_wise(noPos_only, F)
Open_only_pairs <- to_pair_wise(OpenRef_only, F)
DeblurVUnoise_pairs <- to_pair_wise(DeblurVUnoise_format, T)
DadaVUnoise_pairs <- to_pair_wise(DadaVUnoise_format, T)
DadaVOpenref_pairs <- to_pair_wise(DadaVOpenref_format, T)
DeblurVDada_pairs <- to_pair_wise(DeblurVDada_format, T)
DeblurVNoPos_pairs <- to_pair_wise(DeblurVNoPos_format, T)
DeblurVopenref_pairs <- to_pair_wise(DeblurVopenref_format, T)
UnoiseVOpenref_pairs <- to_pair_wise(UnoiseVOpenref_format, T)

#combine all the tables 


#not going to include noPos in the PCoA plot.
#can create another PCoA for that analysis
library(reshape2)
combined_pair_wise <- rbind(Dada_only_pairs, Unoise_only_pairs, Deblur_only_pairs,
                            DeblurVUnoise_pairs, DadaVUnoise_pairs, DadaVOpenref_pairs, DeblurVDada_pairs, 
                            DeblurVopenref_pairs, UnoiseVOpenref_pairs, Open_only_pairs)

combined_pair_wise_cast <- dcast(combined_pair_wise, S1 ~ S2)
rownames(combined_pair_wise_cast) <- combined_pair_wise_cast$S1
combined_pair_wise_cast <- combined_pair_wise_cast[, -1]

for(i in 1:nrow(combined_pair_wise)){
  if(combined_pair_wise[i, "distance"] != combined_pair_wise_cast[as.character(combined_pair_wise[i,"S1"]), 
                                                                  as.character(combined_pair_wise[i,"S2"])]){
    print("failure")
  }
}
head(combined_pair_wise_cast)


#write DM file.
write.table(combined_pair_wise_cast, file="weighted_combined_DM_for_PC.txt", quote = F, sep="\t", col.names=NA)



############################################################unweighted

#load in tables

U_DadaVUnoise <- read.table("DadavUnoise/distances/unweighted_unifrac_rare_DadavUnoise.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = F, check.names=F)
U_DadaVOpenref <- read.table("DadaVOpen-ref/distances/unweighted_unifrac_rare_dadaVopen-ref.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F, check.names=F)
U_DeblurVDada <- read.table("DeblurVDada/distances/unweighted_unifrac_rare_DeblurVDada.txt", sep="\t", header=T, row.names=1, stringsAsFactors = F, check.names=F)
U_DeblurvNoPos <- read.table("DeblurVnoPos/distances/unweighted_unifrac_rare_DeblurVnoPos.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
colnames(U_DeblurvNoPos) <- gsub("001101", "", colnames(U_DeblurvNoPos))
rownames(U_DeblurvNoPos) <- gsub("001101", "", rownames(U_DeblurvNoPos))
U_DeblurVopenref <- read.table("DeblurVOpen-ref/distances/unweighted_unifrac_rare_deblurVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_DeblurVUnoise <- read.table("DeblurVUnoise/distances/unweighted_unifrac_rare_DeblurVUnoise.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_UnoiseVOpenref <- read.table("UnoiseVOpen-ref/distances/unweighted_unifrac_rare_UnoiseVOpen-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_Dada_only <- read.table("Dada/distances/unweighted_unifrac_dada2_header.txt", sep="\t", header=T, row.names = 1, check.names = F, stringsAsFactors = F)
U_Unoise_only <- read.table("Unoise/distances/unweighted_unifrac_unoise_header.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
U_Deblur_only <- read.table("Deblur/distances/unweighted_unifrac_deblur_header.txt", sep="\t", header=T, row.names = 1, check.names = F, stringsAsFactors = F)
U_OpenRef_only <- read.table("Open/distances/unweighted_unifrac_rare_Open-ref.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)
U_noPos_only <- read.table("noPos/distances/unweighted_unifrac_rare_Deblur.txt", sep="\t", header=T, row.names=1, check.names = F, stringsAsFactors = F)

#format the tables
U_DadaVUnoise_format <- format_DM(U_DadaVUnoise, "Dada_", "Unoise_")
U_DadaVOpenref_format <- format_DM(U_DadaVOpenref, "Dada_", "open-ref_")
U_DeblurVDada_format <- format_DM(U_DeblurVDada, "Dada_", "Deblur_")
U_DeblurVNoPos_format <- format_DM(U_DeblurvNoPos, "Deblur_", "Deblur-noPos_")
U_DeblurVopenref_format <- format_DM(U_DeblurVopenref, "Deblur_", "open-ref_")
U_DeblurVUnoise_format <- format_DM(U_DeblurVUnoise, "Deblur_", "Unoise_")
U_UnoiseVOpenref_format <- format_DM(U_UnoiseVOpenref, "Unoise_", "open-ref_")

#to pairwise
U_Dada_only_pairs <- to_pair_wise(U_Dada_only, F)
U_Unoise_only_pairs <- to_pair_wise(U_Unoise_only, F)
U_Deblur_only_pairs <- to_pair_wise(U_Deblur_only, F)
U_noPos_only_pairs <- to_pair_wise(U_noPos_only, F)
U_Open_only_pairs <- to_pair_wise(U_OpenRef_only, F)
U_DeblurVUnoise_pairs <- to_pair_wise(U_DeblurVUnoise_format, T)
U_DadaVUnoise_pairs <- to_pair_wise(U_DadaVUnoise_format, T)
U_DadaVOpenref_pairs <- to_pair_wise(U_DadaVOpenref_format, T)
U_DeblurVDada_pairs <- to_pair_wise(U_DeblurVDada_format, T)
U_DeblurVNoPos_pairs <- to_pair_wise(U_DeblurVNoPos_format, T)
U_DeblurVopenref_pairs <- to_pair_wise(U_DeblurVopenref_format, T)
U_UnoiseVOpenref_pairs <- to_pair_wise(U_UnoiseVOpenref_format, T)

#combined unweighted pairwise tables

#combined the tables
library(reshape2)
U_combined_pair_wise <- rbind(U_Dada_only_pairs, U_Unoise_only_pairs, U_Deblur_only_pairs,
                              U_DeblurVUnoise_pairs, U_DadaVUnoise_pairs, U_DadaVOpenref_pairs, U_DeblurVDada_pairs, 
                              U_DeblurVopenref_pairs, U_UnoiseVOpenref_pairs, U_Open_only_pairs)

U_combined_pair_wise_cast <- dcast(U_combined_pair_wise, S1 ~ S2)
rownames(U_combined_pair_wise_cast) <- U_combined_pair_wise_cast$S1
U_combined_pair_wise_cast <- U_combined_pair_wise_cast[, -1]

for(i in 1:nrow(U_combined_pair_wise)){
  if(U_combined_pair_wise[i, "distance"] != U_combined_pair_wise_cast[as.character(U_combined_pair_wise[i,"S1"]), 
                                                                      as.character(U_combined_pair_wise[i,"S2"])]){
    print("failure")
  }
}


#write DM file.
write.table(U_combined_pair_wise_cast, file="unweighted_combined_DM_for_PC.txt", quote = F, sep="\t", col.names=NA)
