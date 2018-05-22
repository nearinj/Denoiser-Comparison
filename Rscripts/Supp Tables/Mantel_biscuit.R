#Mantel for Biscuit data

library(vegan)


mantel_format <- function(tab, prefix, keep){
  
  ret_tab <- tab[grepl(prefix, rownames(tab)), grepl(prefix, colnames(tab))]
  keep <- lapply(keep, function(x) paste(prefix, x, sep=""))
  ret_tab_fin <- ret_tab[match(keep, rownames(ret_tab)), match(keep, colnames(ret_tab))]
  return(ret_tab_fin)
}

mantel_test <- function(tab, prefix1, prefix2, keep){
  prefix1_tab <- mantel_format(tab, prefix1, keep)
  prefix2_tab <- mantel_format(tab, prefix2, keep)
  prefix1_tab <- prefix1_tab[gsub(prefix2, prefix1, rownames(prefix2_tab)), gsub(prefix2, prefix1, colnames(prefix2_tab))]
  #check if dimensions are the same
  if(!identical(dim(prefix1_tab),dim(prefix2_tab))){return("dims don't match")} 
  
  #if length of true values is equal to total row length then rows are good
  good_to_go_rows <-  length(which(gsub(prefix1, "", rownames(prefix1_tab)) == gsub(prefix2,"", rownames(prefix2_tab)))) == length(rownames(prefix1_tab))
  #if length of true values is equal to total col length then cols are good
  good_to_go_cols <- length(which(gsub(prefix1,"", colnames(prefix1_tab)) == gsub(prefix2,"", colnames(prefix2_tab)))) == length(colnames(prefix2_tab))
  #check if they are the same
  if(good_to_go_rows == T && good_to_go_cols==T){
    return(mantel(prefix1_tab, prefix2_tab, permutations = 999, method="pearson"))
  }else if(good_to_go_rows == F){
    return("Rows don't match")
  }else if(good_to_go_cols == F){
    return("Cols don't match")
  }
}

#read in weighted data


setwd("~/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/fixed_combined/pplacer_distances/")

DM <- read.table("weighted_unifrac_rare_Combined_biscuit_taxa_fixed.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)
rownames(DM) <- gsub("001101", "", rownames(DM))
colnames(DM) <- gsub("001101", "",colnames(DM))

#set up tables to include all of the same samples (Deblur has the smallest number so take the sample list from that)
samples_to_keep <- gsub("Deblur_","",colnames(DM[grepl("Deblur_", rownames(DM)), grepl("Deblur_", colnames(DM))]))


#set mantel format

mantel_DadaVOpenref <- mantel_test(tab = DM, "Dada_", "open-ref_", samples_to_keep)
mantel_DadaVOpenref
mantel_DadaVUnoise <- mantel_test(DM, "Dada_", "Unoise__", samples_to_keep)
mantel_DadaVUnoise
mantel_DeblurVDada <- mantel_test(DM, "Deblur_", "Dada_", samples_to_keep)
mantel_DeblurVDada
mantel_DeblurVnoPos <- mantel_test(DM, "Deblur_", "Deblur-noPos_", samples_to_keep)
mantel_DeblurVnoPos
mantel_DeblurVOpenref <- mantel_test(DM, "Deblur_", "open-ref_", samples_to_keep)
mantel_DeblurVOpenref
mantel_DeblurVUnoise <- mantel_test(DM, "Deblur_", "Unoise__", samples_to_keep)
mantel_DeblurVUnoise
mantel_UnoiseVOpenref <- mantel_test(DM, "Unoise__", "open-ref_", samples_to_keep)
mantel_UnoiseVOpenref


#make r2 table
mantel_r2_table <- data.frame(Dada=double(4),
                              Deblur=double(4),
                              Unoise=double(4),
                              Open_ref_OTU=double(4)
)

rownames(mantel_r2_table) <- c("Dada", "Unoise", "Deblur", "Open_ref_OTU")

mantel_r2_table["Dada", "Unoise"] <- mantel_DadaVUnoise$statistic
mantel_r2_table["Dada", "Deblur"] <- mantel_DeblurVDada$statistic
mantel_r2_table["Dada", "Open_ref_OTU"] <- mantel_DadaVOpenref$statistic
mantel_r2_table["Unoise", "Dada"] <- mantel_DadaVUnoise$statistic
mantel_r2_table["Unoise", "Deblur"] <- mantel_DeblurVUnoise$statistic
mantel_r2_table["Unoise", "Open_ref_OTU"] <- mantel_UnoiseVOpenref$statistic
mantel_r2_table["Deblur", "Dada"] <- mantel_DeblurVDada$statistic
mantel_r2_table["Deblur", "Unoise"] <- mantel_DeblurVUnoise$statistic
mantel_r2_table["Deblur", "Open_ref_OTU"] <- mantel_DeblurVOpenref$statistic
mantel_r2_table["Open_ref_OTU", "Deblur"] <- mantel_DeblurVOpenref$statistic
mantel_r2_table["Open_ref_OTU", "Unoise"] <- mantel_UnoiseVOpenref$statistic
mantel_r2_table["Open_ref_OTU", "Dada"] <-  mantel_DadaVOpenref$statistic
#squre all the values to get R^2 values
mantel_r2_table <- mantel_r2_table[rownames(mantel_r2_table), rownames(mantel_r2_table)]

mantel_p_table <- data.frame(Dada=double(4),
                              Deblur=double(4),
                              Unoise=double(4),
                              Open_ref_OTU=double(4)
)

rownames(mantel_p_table) <- c("Dada", "Unoise", "Deblur", "Open_ref_OTU")


mantel_p_table["Dada", "Unoise"] <- mantel_DadaVUnoise$signif
mantel_p_table["Dada", "Deblur"] <- mantel_DeblurVDada$signif
mantel_p_table["Dada", "Open_ref_OTU"] <- mantel_DadaVOpenref$signif
mantel_p_table["Unoise", "Dada"] <- mantel_DadaVUnoise$signif
mantel_p_table["Unoise", "Deblur"] <- mantel_DeblurVUnoise$signif
mantel_p_table["Unoise", "Open_ref_OTU"] <- mantel_UnoiseVOpenref$signif
mantel_p_table["Deblur", "Dada"] <- mantel_DeblurVDada$signif
mantel_p_table["Deblur", "Unoise"] <- mantel_DeblurVUnoise$signif
mantel_p_table["Deblur", "Open_ref_OTU"] <- mantel_DeblurVOpenref$signif
mantel_p_table["Open_ref_OTU", "Deblur"] <- mantel_DeblurVOpenref$signif
mantel_p_table["Open_ref_OTU", "Unoise"] <- mantel_UnoiseVOpenref$signif
mantel_p_table["Open_ref_OTU", "Dada"] <-  mantel_DadaVOpenref$signif
#squre all the values to get R^2 values
mantel_p_table <- mantel_p_table[rownames(mantel_p_table), rownames(mantel_p_table)]


library(ggplot2)
library(reshape2)
mantel_r2_table[upper.tri(mantel_r2_table)] <- NA


weighted_mantel_melt <- na.omit(melt(mantel_r2_table))

weighted_heat <- ggplot(weighted_mantel_melt, aes(Var1, Var2)) +
  geom_tile(aes(fill=value), color='white') +
  ggtitle(bquote("Mantel "~r^2~"Between Different Weighted Unifrac Distance Matrices in BISCUIT Dataset"))+
  theme_bw() +
  xlab("Method")+ylab("Method")+
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab', limits=c(0,1)) +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank())

weighted_heat


#########################unweighted values

U_DM <- read.table("unweighted_unifrac_rare_Combined_biscuit_taxa_fixed.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors = F)

rownames(U_DM) <- gsub("001101", "", rownames(U_DM))
colnames(U_DM) <- gsub("001101", "", colnames(U_DM))

#set mantel tests
U_mantel_DadaVOpenref <- mantel_test(tab = U_DM, "Dada_", "open-ref_", samples_to_keep)
U_mantel_DadaVOpenref
U_mantel_DadaVUnoise <- mantel_test(U_DM, "Dada_", "Unoise__", samples_to_keep)
U_mantel_DadaVUnoise
U_mantel_DeblurVDada <- mantel_test(U_DM, "Deblur_", "Dada_", samples_to_keep)
U_mantel_DeblurVDada
U_mantel_DeblurVnoPos <- mantel_test(U_DM, "Deblur_", "Deblur-noPos_", samples_to_keep)
U_mantel_DeblurVnoPos
U_mantel_DeblurVOpenref <- mantel_test(U_DM, "Deblur_", "open-ref_", samples_to_keep)
U_mantel_DeblurVOpenref
U_mantel_DeblurVUnoise <- mantel_test(U_DM, "Deblur_", "Unoise__", samples_to_keep)
U_mantel_DeblurVUnoise
U_mantel_UnoiseVOpenref <- mantel_test(U_DM, "Unoise__", "open-ref_", samples_to_keep)
U_mantel_UnoiseVOpenref

U_mantel_r2_table <- data.frame(Dada=double(4),
                              Deblur=double(4),
                              Unoise=double(4),
                              Open_ref_OTU=double(4)
)

rownames(U_mantel_r2_table) <- c("Dada", "Unoise", "Deblur", "Open_ref_OTU")

U_mantel_r2_table["Dada", "Unoise"] <- U_mantel_DadaVUnoise$statistic
U_mantel_r2_table["Dada", "Deblur"] <- U_mantel_DeblurVDada$statistic
U_mantel_r2_table["Dada", "Open_ref_OTU"] <- U_mantel_DadaVOpenref$statistic
U_mantel_r2_table["Unoise", "Dada"] <- U_mantel_DadaVUnoise$statistic
U_mantel_r2_table["Unoise", "Deblur"] <- U_mantel_DeblurVUnoise$statistic
U_mantel_r2_table["Unoise", "Open_ref_OTU"] <- U_mantel_UnoiseVOpenref$statistic
U_mantel_r2_table["Deblur", "Dada"] <- U_mantel_DeblurVDada$statistic
U_mantel_r2_table["Deblur", "Unoise"] <- U_mantel_DeblurVUnoise$statistic
U_mantel_r2_table["Deblur", "Open_ref_OTU"] <- U_mantel_DeblurVOpenref$statistic
U_mantel_r2_table["Open_ref_OTU", "Deblur"] <- U_mantel_DeblurVOpenref$statistic
U_mantel_r2_table["Open_ref_OTU", "Unoise"] <- U_mantel_UnoiseVOpenref$statistic
U_mantel_r2_table["Open_ref_OTU", "Dada"] <-  U_mantel_DadaVOpenref$statistic
#squre all the values to get R^2 values
U_mantel_r2_table <- U_mantel_r2_table[rownames(U_mantel_r2_table), rownames(U_mantel_r2_table)]



U_mantel_p_table <- data.frame(Dada=double(4),
                                Deblur=double(4),
                                Unoise=double(4),
                                Open_ref_OTU=double(4)
)

rownames(U_mantel_p_table) <- c("Dada", "Unoise", "Deblur", "Open_ref_OTU")

U_mantel_p_table["Dada", "Unoise"] <- U_mantel_DadaVUnoise$signif
U_mantel_p_table["Dada", "Deblur"] <- U_mantel_DeblurVDada$signif
U_mantel_p_table["Dada", "Open_ref_OTU"] <- U_mantel_DadaVOpenref$signif
U_mantel_p_table["Unoise", "Dada"] <- U_mantel_DadaVUnoise$signif
U_mantel_p_table["Unoise", "Deblur"] <- U_mantel_DeblurVUnoise$signif
U_mantel_p_table["Unoise", "Open_ref_OTU"] <- U_mantel_UnoiseVOpenref$signif
U_mantel_p_table["Deblur", "Dada"] <- U_mantel_DeblurVDada$signif
U_mantel_p_table["Deblur", "Unoise"] <- U_mantel_DeblurVUnoise$signif
U_mantel_p_table["Deblur", "Open_ref_OTU"] <- U_mantel_DeblurVOpenref$signif
U_mantel_p_table["Open_ref_OTU", "Deblur"] <- U_mantel_DeblurVOpenref$signif
U_mantel_p_table["Open_ref_OTU", "Unoise"] <- U_mantel_UnoiseVOpenref$signif
U_mantel_p_table["Open_ref_OTU", "Dada"] <-  U_mantel_DadaVOpenref$signif
#squre all the values to get R^2 values
U_mantel_p_table <- U_mantel_p_table[rownames(U_mantel_p_table), rownames(U_mantel_p_table)]


U_mantel_r2_table[upper.tri(U_mantel_r2_table)] <- NA


unweighted_mantel_melt <- na.omit(melt(U_mantel_r2_table))

unweighted_heat <- ggplot(unweighted_mantel_melt, aes(Var1, Var2)) +
  geom_tile(aes(fill=value), color='white') +
  ggtitle(bquote("Mantel "~r^2~"Between Different Unweighted Unifrac Distance Matrices in BISCUIT Dataset"))+
  theme_bw() +
  xlab("Method")+ylab("Method")+
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab', limits=c(0,1)) +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank())

unweighted_heat


##################################BRAY CURTIS#####################################################


B_DM <- read.table("../../final_combined/bray_curtis/distances/bray_curtis_rare_CombinedV2_biscuit_taxa_L6.txt", sep="\t", row.names=1, check.names=F, stringsAsFactors = F, header=T)


rownames(B_DM) <- gsub("001101", "", rownames(B_DM))
colnames(B_DM) <- gsub("001101", "", colnames(B_DM))

#set mantel tests
B_mantel_DadaVOpenref <- mantel_test(tab = B_DM, "Dada_", "open-ref_", samples_to_keep)
B_mantel_DadaVOpenref
B_mantel_DadaVUnoise <- mantel_test(B_DM, "Dada_", "Unoise__", samples_to_keep)
B_mantel_DadaVUnoise
B_mantel_DeblurVDada <- mantel_test(B_DM, "Deblur_", "Dada_", samples_to_keep)
B_mantel_DeblurVDada
B_mantel_DeblurVnoPos <- mantel_test(B_DM, "Deblur_", "Deblur-noPos_", samples_to_keep)
B_mantel_DeblurVnoPos
B_mantel_DeblurVOpenref <- mantel_test(B_DM, "Deblur_", "open-ref_", samples_to_keep)
B_mantel_DeblurVOpenref
B_mantel_DeblurVUnoise <- mantel_test(B_DM, "Deblur_", "Unoise__", samples_to_keep)
B_mantel_DeblurVUnoise
B_mantel_UnoiseVOpenref <- mantel_test(B_DM, "Unoise__", "open-ref_", samples_to_keep)
B_mantel_UnoiseVOpenref

B_mantel_r2_table <- data.frame(Dada=double(4),
                                Deblur=double(4),
                                Unoise=double(4),
                                Open_ref_OTU=double(4)
)

rownames(B_mantel_r2_table) <- c("Dada", "Unoise", "Deblur", "Open_ref_OTU")

B_mantel_r2_table["Dada", "Unoise"] <- B_mantel_DadaVUnoise$statistic
B_mantel_r2_table["Dada", "Deblur"] <- B_mantel_DeblurVDada$statistic
B_mantel_r2_table["Dada", "Open_ref_OTU"] <- B_mantel_DadaVOpenref$statistic
B_mantel_r2_table["Unoise", "Dada"] <- B_mantel_DadaVUnoise$statistic
B_mantel_r2_table["Unoise", "Deblur"] <- B_mantel_DeblurVUnoise$statistic
B_mantel_r2_table["Unoise", "Open_ref_OTU"] <- B_mantel_UnoiseVOpenref$statistic
B_mantel_r2_table["Deblur", "Dada"] <- B_mantel_DeblurVDada$statistic
B_mantel_r2_table["Deblur", "Unoise"] <- B_mantel_DeblurVUnoise$statistic
B_mantel_r2_table["Deblur", "Open_ref_OTU"] <- B_mantel_DeblurVOpenref$statistic
B_mantel_r2_table["Open_ref_OTU", "Deblur"] <- B_mantel_DeblurVOpenref$statistic
B_mantel_r2_table["Open_ref_OTU", "Unoise"] <- B_mantel_UnoiseVOpenref$statistic
B_mantel_r2_table["Open_ref_OTU", "Dada"] <-  B_mantel_DadaVOpenref$statistic
#squre all the values to get R^2 values
B_mantel_r2_table <- B_mantel_r2_table[rownames(B_mantel_r2_table), rownames(B_mantel_r2_table)]^2
