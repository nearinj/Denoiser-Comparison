#Script to make Figure 4

#first load in total sequences for each dataset
library(ggplot2)

Blueberry_sequences <- c(11613, 8270, 16609, 21297)
Biscuit_sequences <- c(1160, 1086, 1589, 1898)
Exercise_sequences <- c(1643, 1318, 1927, 6336)

Methods <- c("UNOISE3", "Deblur", "DADA2", "OTU")

Blueberry <- as.data.frame(cbind(Blueberry_sequences, Methods))
Blueberry$Methods <- factor(Blueberry$Methods, levels=Blueberry$Methods)
Blueberry$Blueberry_sequences <- as.numeric(as.character(Blueberry$Blueberry_sequences))

Biscuit <- as.data.frame(cbind(Biscuit_sequences, Methods))
Biscuit$Methods <- factor(Biscuit$Methods, levels=Biscuit$Methods)
Biscuit$Biscuit_sequences <- as.numeric(as.character(Biscuit$Biscuit_sequences))

Exercise <- as.data.frame(cbind(Exercise_sequences, Methods))
Exercise$Methods <- factor(Exercise$Methods, levels=Exercise$Methods)
Exercise$Exercise_sequences <- as.numeric(as.character(Exercise$Exercise_sequences))

#colors for the groupings
group.colors <- c(Deblur = "#333BFF", DADA2 = "#CC6600", UNOISE3 = "#9633FF", OTU = "#ff6666")

#plot the graphs
Blueberry_ASV_count <- ggplot(Blueberry, aes(x=Methods, y=Blueberry_sequences, fill=Methods)) + geom_bar(stat="identity") + ylab("Total ASVs/OTUs") + 
  scale_fill_manual(values=group.colors)
Blueberry_ASV_count

Biscuit_ASV_count <- ggplot(Biscuit, aes(x=Methods, y=Biscuit_sequences, fill=Methods)) + geom_bar(stat='identity') + ylab("Total ASVs/OTUs")+ 
  scale_fill_manual(values=group.colors)
Biscuit_ASV_count

Exercise_ASV_count <- ggplot(Exercise, aes(x=Methods, y=Exercise_sequences, fill=Methods)) + geom_bar(stat='identity') + ylab("Total ASVs/OTUs")+
  scale_fill_manual(values=group.colors)
Exercise_ASV_count


#make heat maps from correlations of otu numbers
#functions
#generates the vals column in the heatmap table with Rho and * if significant
get_significance <- function(x){
  ret_list <- list()
  j <- 1
  for ( i in 1:nrow(x)){
    
    lab <- sprintf("%0.3f", x$Rho[i])
    
    if(x$p[i] < 0.0005){
      ret_list[j] <- paste(lab, "*", sep=" ")
    }else if(x$p[i] < 0.005){
      ret_list[j] <- paste(lab, "*", sep=" ")
    }else if(x$p[i] < 0.05){
      ret_list[j] <- paste(lab, "*", sep=" ")
    }else{
      ret_list[j] <- paste(lab)
    }
    
    j <- j +1
  }
  return(ret_list)
}

#generates the table that the heatmap will be built out of

generate_table_heatmap <- function(x) {
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
  blue_deblur_dada <- cor.test(blue_dada$observed_otus, blue_deblur$observed_otus, method='spearman')
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

setwd("~/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/fixed_combined/")

alpha_data_blue <- read.table("rare_alpha_div.tsv", sep="\t", header=T) 
blue_heat_table <- generate_table_heatmap(alpha_data_blue)



blue_heat <- ggplot(blue_heat_table, aes(M1, M2, fill=Rho)) + geom_tile() + scale_fill_gradient(low="yellow", high="red", limits=c(0,1)) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  geom_text(data=blue_heat_table, aes(label=vals))
  


blue_heat



#heatmap for biscuit data
setwd("~/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/fixed_combined/")

biscuit_alpha <- read.table("rare_alpha_div.tsv", header=T, sep="\t")

biscuit_heat_table <- generate_table_heatmap(biscuit_alpha)


biscuit_heat <- ggplot(biscuit_heat_table, aes(M1, M2, fill=Rho)) + geom_tile() + scale_fill_gradient(low="yellow", high="red", limits=c(0,1)) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  geom_text(data=biscuit_heat_table, aes(label=vals))



biscuit_heat


#heatmap for Exercise data
setwd("~/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/fixed_combined/")

exercise_alpha <- read.table("rare_alpha_div.tsv", sep="\t", header=T)

exercise_heat_table <- generate_table_heatmap(exercise_alpha)


exercise_heat <- ggplot(exercise_heat_table, aes(M1, M2, fill=Rho)) + geom_tile() + scale_fill_gradient(low="yellow", high="red", limits=c(0,1)) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  geom_text(data=exercise_heat_table, aes(label=vals))



exercise_heat


#need to do some bug checking first!
#combined all the plots together
library(cowplot)

fig4 <- plot_grid(Blueberry_ASV_count, blue_heat, Exercise_ASV_count, exercise_heat, Biscuit_ASV_count, biscuit_heat, labels="AUTO", nrow=3)
fig4


#make a test dataframe to test the heatmap function





