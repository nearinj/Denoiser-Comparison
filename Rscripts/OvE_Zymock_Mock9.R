library(ggplot2)
library(cowplot)

setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")


Zymock_DADA_HIGH <- makeFullTable("per_Zymock_Dada_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DADA_HIGH <- CalculateOvE_Zymock(Zymock_DADA_HIGH)
Zymock_DEBLUR_HIGH <- makeFullTable("per_Zymock_Deblur_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_HIGH <- CalculateOvE_Zymock(Zymock_DEBLUR_HIGH)
Zymock_UNOISE_HIGH <- makeFullTable("per_Zymock_Unoise_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_HIGH <- CalculateOvE_Zymock(Zymock_UNOISE_HIGH)

GraphOvE_Zymock(Zymock_DADA_HIGH, Zymock_DEBLUR_HIGH, Zymock_UNOISE_HIGH, "Zymock_High")

Zymock_DADA_MED <- makeFullTable("per_Zymock_Dada_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DADA_MED <- CalculateOvE_Zymock(Zymock_DADA_MED)
Zymock_DEBLUR_MED <- makeFullTable("per_Zymock_Deblur_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_MED <- CalculateOvE_Zymock(Zymock_DEBLUR_MED)
Zymock_UNOISE_MED <- makeFullTable("per_Zymock_Unoise_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_MED <- CalculateOvE_Zymock(Zymock_UNOISE_MED)

GraphOvE_Zymock(Zymock_DADA_MED, Zymock_DEBLUR_MED, Zymock_UNOISE_MED, "Zymock_Med")

Zymock_DADA_LOW <- makeFullTable("per_Zymock_Dada_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DADA_LOW <- CalculateOvE_Zymock(Zymock_DADA_LOW)
Zymock_DEBLUR_LOW <- makeFullTable("per_Zymock_Deblur_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_LOW <- CalculateOvE_Zymock(Zymock_DEBLUR_LOW)
Zymock_UNOISE_LOW <- makeFullTable("per_Zymock_Unoise_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_LOW <- CalculateOvE_Zymock(Zymock_UNOISE_LOW)

GraphOvE_Zymock(Zymock_DADA_LOW, Zymock_DEBLUR_LOW, Zymock_UNOISE_LOW, "Zymock_Low")

Mock9_DADA_HIGH <- makeFullTable("per_Mock-9_Dada_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_HIGH <- CalculateOvE_Zymock(Mock9_DADA_HIGH)
Mock9_DEBLUR_HIGH <- makeFullTable("per_Mock-9_Deblur_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_HIGH <- CalculateOvE_Zymock(Mock9_DEBLUR_HIGH)
Mock9_UNOISE_HIGH <- makeFullTable("per_Mock-9_Unoise_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_HIGH <- CalculateOvE_Zymock(Mock9_UNOISE_HIGH)

GraphOvE_Zymock(Mock9_DADA_HIGH, Mock9_DEBLUR_HIGH, Mock9_UNOISE_HIGH, "Mock-9_High")

Mock9_DADA_MED <- makeFullTable("per_Mock-9_Dada_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_MED <- CalculateOvE_Zymock(Mock9_DADA_MED)
Mock9_DEBLUR_MED <- makeFullTable("per_Mock-9_Deblur_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_MED <- CalculateOvE_Zymock(Mock9_DEBLUR_MED)
Mock9_UNOISE_MED <- makeFullTable("per_Mock-9_Unoise_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_MED <- CalculateOvE_Zymock(Mock9_UNOISE_MED)

GraphOvE_Zymock(Mock9_DADA_MED, Mock9_DEBLUR_MED, Mock9_UNOISE_MED, "Mock-9_Med")

Mock9_DADA_LOW <- makeFullTable("per_Mock-9_Dada_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_LOW <- CalculateOvE_Zymock(Mock9_DADA_LOW)
Mock9_DEBLUR_LOW <- makeFullTable("per_Mock-9_Deblur_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_LOW <- CalculateOvE_Zymock(Mock9_DEBLUR_LOW)
Mock9_UNOISE_LOW <- makeFullTable("per_Mock-9_Unoise_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_LOW <- CalculateOvE_Zymock(Mock9_UNOISE_LOW)

GraphOvE_Zymock(Mock9_DADA_LOW, Mock9_DEBLUR_LOW, Mock9_UNOISE_LOW, "Mock-9_Low")

makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  return(final)
}

CalculateOvE_Zymock <- function(tab){
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
  
  
  vars <- c("Sample1OvE", "Sample2OvE", "Sample3OvE")
  temptab <- tab[vars]
  means <- rowMeans(temptab)
  tab$OvEMeans <- means
  
  return(tab)
}

GraphOvE_Zymock <- function(tab1,tab2,tab3,filt){
  
  #table 1
  #remove total row
  tab1_filt <- tab1[-grep("Total", tab1$Organism),]
  #melt data into long format to use in ggplot
  tab1_melt <- melt(data=tab1_filt, measure.vars=c("Sample1OvE", "Sample2OvE", "Sample3OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plot1 <- ggplot(tab1_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.55))
  
  # table 2
  #remove total row
  tab2_filt <- tab2[-grep("Total", tab2$Organism),]
  #melt data into long format to use in ggplot
  tab2_melt <- melt(data=tab2_filt, measure.vars=c("Sample1OvE", "Sample2OvE", "Sample3OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  plot2 <- ggplot(tab2_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.55))
  
  #table 3
  #remove total row
  tab3_filt <- tab3[-grep("Total", tab3$Organism),]
  #melt data into long format to use in ggplot
  tab3_melt <- melt(data=tab3_filt, measure.vars=c("Sample1OvE", "Sample2OvE", "Sample3OvE"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plot3 <- ggplot(tab3_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.55))
  
  
  
  
  #table 4 the means of each sample!
  names(tab1)[names(tab1) == 'OvEMeans'] <- 'Dada'
  names(tab2)[names(tab2) == 'OvEMeans'] <- 'Deblur'
  names(tab3)[names(tab3) == 'OvEMeans'] <- 'Unoise'
  
  
  
  temptabM <- data.frame(Organism=tab1$Organism, Dada=tab1$Dada, Deblur=tab2$Deblur, Unoise=tab3$Unoise)
  
  
  tabM_filt <- temptabM[-grep("Total", temptabM$Organism),]
  #melt data into long format to use in ggplot
  tabM_melt <- melt(data=tabM_filt, measure.vars=c("Dada", "Deblur", "Unoise"), 
                    id.vars = c("Organism"), value.name="OvE", variable.name="sample")
  
  plotM <- ggplot(tabM_melt, aes(x=sample, y=OvE)) +theme(axis.text.x=element_text(size=10)) + geom_boxplot(outlier.colour = "red", outlier.shape = 8,) + scale_y_continuous(limits = c(0,.55))
  
  
  
  grid_plot <- plot_grid(plot1, plot2, plot3, plotM, labels=c("A", "B", "C", "D"))
  filt <- paste("../OvE", filt, sep="/")
  save_plot(paste(filt, "png", sep ="."), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 1.3)
  
  
  
}
