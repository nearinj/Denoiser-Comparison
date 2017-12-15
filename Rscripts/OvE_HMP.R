#Used to make the Observed vs Expected for the Abundances of the HMP study
library(ggplot2)
library(cowplot)

setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")

HMP_DADA_HIGH <- makeFullTable("per_HMP_Dada_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_HIGH <- CalculateOvE_HMP(HMP_DADA_HIGH)

HMP_DEBLUR_HIGH <- makeFullTable("per_HMP_Deblur_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_HIGH <- CalculateOvE_HMP(HMP_DEBLUR_HIGH)

HMP_UNOISE_HIGH <- makeFullTable("per_HMP_Unoise_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_HIGH <- CalculateOvE_HMP(HMP_UNOISE_HIGH)

GraphOvE_HMP(HMP_DADA_HIGH, HMP_DEBLUR_HIGH, HMP_UNOISE_HIGH, "HMP_High")

HMP_DADA_MED <- makeFullTable("per_HMP_Dada_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_MED <- CalculateOvE_HMP(HMP_DADA_MED)
HMP_DEBLUR_MED <- makeFullTable("per_HMP_Deblur_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_MED <- CalculateOvE_HMP(HMP_DEBLUR_MED)
HMP_UNOISE_MED <- makeFullTable("per_HMP_Unoise_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_MED <- CalculateOvE_HMP(HMP_UNOISE_MED)

GraphOvE_HMP(HMP_DADA_MED, HMP_DEBLUR_MED, HMP_UNOISE_MED, "HMP_Med")

HMP_DADA_LOW <- makeFullTable("per_HMP_Dada_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_LOW <- CalculateOvE_HMP(HMP_DADA_LOW)
HMP_DEBLUR_LOW <- makeFullTable("per_HMP_Deblur_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_LOW <- CalculateOvE_HMP(HMP_DEBLUR_LOW)
HMP_UNOISE_LOW <- makeFullTable("per_HMP_Unoise_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_LOW <- CalculateOvE_HMP(HMP_UNOISE_LOW)

GraphOvE_HMP(HMP_DADA_LOW, HMP_DEBLUR_LOW, HMP_UNOISE_LOW, "HMP_Low")


makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  return(final)
}


CalculateOvE_HMP <- function(tab){
  #fill in any Na's with 0s in the samples!
  tab[is.na(tab)] <- 0
  
  #get OvE
  tab$Sample1OvE <- (tab$Sample1-tab$Expected)
  tab$Sample1OvE <- tab$Sample1OvE^2
  tab$Sample1OvE <- tab$Sample1OvE^0.5
  
  tab$Sample2OvE <- (tab$Sample2-tab$Expected)
  tab$Sample2OvE <- tab$Sample2OvE^2
  tab$Sample2OvE <- tab$Sample2OvE^0.5
  
  tab$Sample3OvE <- (tab$Sample3-tab$Expected)
  tab$Sample3OvE <- tab$Sample3OvE^2
  tab$Sample3OvE <- tab$Sample3OvE^0.5
  
  
  tab$Sample4OvE <- (tab$Sample4-tab$Expected)
  tab$Sample4OvE <- tab$Sample4OvE^2
  tab$Sample4OvE <- tab$Sample4OvE^0.5

  vars <- c("Sample1OvE", "Sample2OvE", "Sample3OvE", "Sample4OvE")
  temptab <- tab[vars]
  means <- rowMeans(temptab)
  tab$OvEMeans <- means
  
  return(tab)
}




GraphOvE_HMP <- function(tab1,tab2,tab3,filt){
  
  #table 1
   #remove total row
  tab1_filt <- tab1[-grep("Total", tab1$Organism),]
  #melt data into long format to use in ggplot
  tab1_melt <- melt(data=tab1_filt, measure.vars=c("Sample1OvE", "Sample2OvE", "Sample3OvE", "Sample4OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plot1 <- ggplot(tab1_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits=c(0,.1))
  
  # table 2
  #remove total row
  tab2_filt <- tab2[-grep("Total", tab2$Organism),]
  #melt data into long format to use in ggplot
  tab2_melt <- melt(data=tab2_filt, measure.vars=c("Sample1OvE", "Sample2OvE", "Sample3OvE", "Sample4OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  plot2 <- ggplot(tab2_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits=c(0,.1))
  
  #table 3
  #remove total row
  tab3_filt <- tab3[-grep("Total", tab3$Organism),]
  #melt data into long format to use in ggplot
  tab3_melt <- melt(data=tab3_filt, measure.vars=c("Sample1OvE", "Sample2OvE", "Sample3OvE", "Sample4OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plot3 <- ggplot(tab3_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits=c(0,.1))

  
  
  
  #table 4 the means of each sample!
  names(tab1)[names(tab1) == 'OvEMeans'] <- 'Dada'
  names(tab2)[names(tab2) == 'OvEMeans'] <- 'Deblur'
  names(tab3)[names(tab3) == 'OvEMeans'] <- 'Unoise'
  
  
  temptabM <- data.frame(Organism=tab1$Organism, Dada=tab1$Dada, Deblur=tab2$Deblur, Unoise=tab3$Unoise)
  
  
  tabM_filt <- temptabM[-grep("Total", temptabM$Organism),]
  #melt data into long format to use in ggplot
  tabM_melt <- melt(data=tabM_filt, measure.vars=c("Dada", "Deblur", "Unoise"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plotM <- ggplot(tabM_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits=c(0,.1))
  
  
  
  grid_plot <- plot_grid(plot1, plot2, plot3, plotM, labels=c("A", "B", "C", "D"))
  filt <- paste("../OvE", filt, sep="/")
  save_plot(paste(filt, "png", sep ="."), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 1.3)
  
  
  
}




