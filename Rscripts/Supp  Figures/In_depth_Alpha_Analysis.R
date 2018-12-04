#Further analysis of alpha data



setwd("~/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/fixed_combined/")



#read in alpha data

alpha_data_blue <- read.table("testing/test.tsv", sep="\t", header=T) 

#need to edit funcation so u can input name of column u want to gen heatmap for.

generate_table_heatmap <- function(x, colname) {
  colnames(x)[1] <- "Method"
  x$sample <- gsub(".*_","",x$Method)
  x$Method <- gsub("_.*","",x$Method)
  
  
  #split them and order them all the same to do correlation test
  blue_deblur <- x[which(x$Method=="Deblur"),]
  rownames(blue_deblur) <- blue_deblur$sample
  blue_dada <- x[which(x$Method=="Dada"),]
  rownames(blue_dada) <- blue_dada$sample
  blue_unoise <- x[which(x$Method=="Unoise"),]
  rownames(blue_unoise) <- blue_unoise$sample
  blue_OTU <- x[which(x$Method=="open-ref"),]
  rownames(blue_OTU) <- blue_OTU$sample
  
  blue_deblur <- blue_deblur[rownames(blue_dada),]
  blue_unoise <- blue_unoise[rownames(blue_dada),]
  blue_OTU <- blue_OTU[rownames(blue_dada),]
  
  #check to make sure that they are all in the same order.
  print(length(which(blue_deblur$sample == blue_dada$sample)))
  print(length(which(blue_unoise$sample == blue_dada$sample)))
  print(length(which(blue_OTU$sample == blue_dada$sample)))
  
  
  
  #do the correlaton tests
  #deblurVdada
  blue_deblur_dada <- cor.test('$'(blue_dada, colname), '$'(blue_deblur, colname), method='spearman')
  blue_deblur_dada_stats <- c("DADA2", "Deblur", blue_deblur_dada$estimate, blue_deblur_dada$p.value)
  #DeblurVUnoise
  blue_deblur_unoise <- cor.test(blue_deblur$observed_otus, blue_unoise$observed_otus, method='spearman')
  blue_deblur_unoise_stats <- c("UNOISE3", "Deblur", blue_deblur_unoise$estimate, blue_deblur_unoise$p.value)
  #DadaVUnoise
  blue_dada_unoise <- cor.test(blue_dada$observed_otus, blue_unoise$observed_otus, method='spearman')
  blue_dada_unoise_stats <- c("DADA2", "UNOISE3", blue_dada_unoise$estimate, blue_dada_unoise$p.value)
  #OtuVdada
  blue_dada_OTU <- cor.test(blue_dada$observed_otus, blue_OTU$observed_otus, method='spearman')
  blue_dada_OTU_stats <- c("DADA2", "OTU", blue_dada_OTU$estimate, blue_dada_OTU$p.value)
  #OTUVDeblur
  blue_deblur_OTU <- cor.test(blue_deblur$observed_otus, blue_OTU$observed_otus, method='spearman')
  blue_deblur_OTU_stats <- c("Deblur", "OTU", blue_deblur_OTU$estimate, blue_deblur_OTU$p.value)
  #OTUVUnnoise
  blue_unoise_OTU <- cor.test(blue_unoise$observed_otus, blue_OTU$observed_otus, method='spearman')
  blue_unoise_OTU_stats <- c("UNOISE3", "OTU", blue_unoise_OTU$estimate, blue_unoise_OTU$p.value)
  
  
  #start making heat map
  #hmmmmm make 
  #R value table
  #make a melted table to contain all the values
  #i have all the test statistics I want to make a heatmap based on R value and then put * by those that are significant...
  #y axis order
  Methodsy<- c("OTU", "Deblur", "UNOISE3", "DADA2")
  #x axis order
  Methodsx <- c("DADA2", "UNOISE3", "Deblur", "OTU")
  #build table
  blue_cors <- as.data.frame(rbind(blue_deblur_dada_stats, #Deblur V DADA
                                   blue_deblur_unoise_stats, #Deblur V UNOISE
                                   blue_dada_unoise_stats, #Dada V UNOISE
                                   blue_dada_OTU_stats, #DAda V OTU
                                   blue_deblur_OTU_stats, #Deblur V OTU
                                   blue_unoise_OTU_stats)) #unoise V OTU
  colnames(blue_cors) <- c("M1", "M2", "Rho", "p")
  blue_cors$M2 <- factor(blue_cors$M2, levels=Methodsy)
  blue_cors$M1 <- factor(blue_cors$M1, levels=Methodsx)
  blue_cors$Rho <- as.numeric(as.character(blue_cors$Rho))
  blue_cors$p <- as.numeric(as.character(blue_cors$p))
  
  #get vals with Rho value + * if sig
  blue_cors$vals <- as.character(get_significance(blue_cors))
  
  return(blue_cors)
}

test <- generate_table_heatmap(alpha_data_blue)
