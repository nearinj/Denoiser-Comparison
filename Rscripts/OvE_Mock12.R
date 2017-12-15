#Used to make the observed vs expected for the Mock12 study
library(ggplot2)
library(cowplot)
setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")

Mock12_DADA_HIGH <- makeFullTable("per_Mock-12_Dada_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_HIGH <- CalculateOvE_Mock12(Mock12_DADA_HIGH)
Mock12_DEBLUR_HIGH <- makeFullTable("per_Mock-12_Deblur_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_HIGH <- CalculateOvE_Mock12(Mock12_DEBLUR_HIGH)
Mock12_UNOISE_HIGH <- makeFullTable("per_Mock-12_Unoise_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOISE_HIGH <- CalculateOvE_Mock12(Mock12_UNOISE_HIGH)

GraphOvE_Mock12(Mock12_DADA_HIGH, Mock12_DEBLUR_HIGH, Mock12_UNOISE_HIGH, "Mock-12_High")

Mock12_DADA_MED <- makeFullTable("per_Mock-12_Dada_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_MED <- CalculateOvE_Mock12(Mock12_DADA_MED)
Mock12_DEBLUR_MED <- makeFullTable("per_Mock-12_Deblur_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_MED <- CalculateOvE_Mock12(Mock12_DEBLUR_MED)
Mock12_UNOISE_MED <- makeFullTable("per_Mock-12_Unoise_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOISE_MED <- CalculateOvE_Mock12(Mock12_UNOISE_MED)

GraphOvE_Mock12(Mock12_DADA_MED, Mock12_DEBLUR_MED, Mock12_UNOISE_MED, "Mock-12_Med")

Mock12_DADA_LOW <- makeFullTable("per_Mock-12_Dada_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_LOW <- CalculateOvE_Mock12(Mock12_DADA_LOW)
Mock12_DEBLUR_LOW <- makeFullTable("per_Mock-12_Deblur_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_LOW <- CalculateOvE_Mock12(Mock12_DEBLUR_LOW)
Mock12_UNOISE_LOW <- makeFullTable("per_Mock-12_Unoise_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOISE_LOW <- CalculateOvE_Mock12(Mock12_UNOISE_LOW)

GraphOvE_Mock12(Mock12_DADA_LOW, Mock12_DEBLUR_LOW, Mock12_UNOISE_LOW, "Mock-12_Low")


makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  return(final)
}

CalculateOvE_Mock12 <- function(tab){
  #fill in any Na's with 0s in the samples!
  tab[is.na(tab)] <- 0
  
  #get OvE
  tab$Sample1OvE <- (tab$Sample1-tab$Expected)
  tab$Sample1OvE <- tab$Sample1OvE^2
  tab$Sample1OvE <- tab$Sample1OvE^0.5

  vars <- c("Sample1OvE")
  temptab <- tab[vars]
  means <- rowMeans(temptab)
  tab$OvEMeans <- means
  
  return(tab)
}


GraphOvE_Mock12 <- function(tab1,tab2,tab3,filt){
  
  #table 1
  #remove total row
  tab1_filt <- tab1[-grep("Total", tab1$Organism),]
  #melt data into long format to use in ggplot
  tab1_melt <- melt(data=tab1_filt, measure.vars=c("Sample1OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plot1 <- ggplot(tab1_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.1))
  
  # table 2
  #remove total row
  tab2_filt <- tab2[-grep("Total", tab2$Organism),]
  #melt data into long format to use in ggplot
  tab2_melt <- melt(data=tab2_filt, measure.vars=c("Sample1OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  plot2 <- ggplot(tab2_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.1))
  
  #table 3
  #remove total row
  tab3_filt <- tab3[-grep("Total", tab3$Organism),]
  #melt data into long format to use in ggplot
  tab3_melt <- melt(data=tab3_filt, measure.vars=c("Sample1OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plot3 <- ggplot(tab3_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.1))
  
  
  
  
  #table 4 the means of each sample!
  names(tab1)[names(tab1) == 'OvEMeans'] <- 'Dada'
  names(tab2)[names(tab2) == 'OvEMeans'] <- 'Deblur'
  names(tab3)[names(tab3) == 'OvEMeans'] <- 'Unoise'
  
  
  
  temptabM <- data.frame(Organism=tab1$Organism, Dada=tab1$Dada, Deblur=tab2$Deblur, Unoise=tab3$Unoise)
  
  
  tabM_filt <- temptabM[-grep("Total", temptabM$Organism),]
  #melt data into long format to use in ggplot
  tabM_melt <- melt(data=tabM_filt, measure.vars=c("Dada", "Deblur", "Unoise"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plotM <- ggplot(tabM_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.1))
  
  
  
  grid_plot <- plot_grid(plot1, plot2, plot3, plotM, labels=c("A", "B", "C", "D"))
  filt <- paste("../OvE", filt, sep="/")
  save_plot(paste(filt, "png", sep ="."), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 1.3)
  
  
  
}
