### similarity propotion curve for mock communities

library(ggplot2)
library(cowplot)
library(gridGraphics)
setwd("~/projects/DenoiseCompare/Data_Analysis/fixed_ASV_analysis/Sims/")
#name list for the x lab for each graph
Pipeline_name_list <- c("DADA2", "Deblur", "UNOISE3", "Open")

#function to make histograms

plot_prop_curve <- function(list, namelist, ylim){
  plot_list=list()
  for(i in 1:length(list)){
    #set anything below 75 to be 75 for easy binning.... is this a great cut off point? 
    list[[i]]$V2[list[[i]]$V2 < 75] <- 71
    print(length(which(list[[i]]$V2 == 71)))
    plot_list[[i]] <- ggplot(list[[i]], aes(V2)) + geom_histogram(breaks=c(seq(70,100, by=1))) + ylim(c(0,ylim[i]))+ xlab(namelist[[i]])+
      scale_x_continuous(limits=c(70,100), breaks=c(seq(70, 100, by=5)), labels = c("< 75", seq(75,100,by=5)))
  }
  res <- plot_grid(plotlist=plot_list, labels="AUTO")
}


#High filt data

#read in only everyone second line (first line is the sequence name, which is not needed for this plot)

#Dumby test data
V1 <- c(1,2,3,4,5,6,7,8,9,10)

V2 <- c(0,97,90,74,90,100,99,99.9,86,75)
testdf <- as.data.frame(V1)
testdf$V2 <- V2
testlist <- list(testdf, testdf)
test_plot <- plot_prop_curve(testlist, Pipeline_name_list, ylim=100)
test_plot

#HMP
Dada_high_HMP <- read.table("HMP_DadaHigh.tsv", header=F, sep="\t")
Deblur_high_HMP <- read.table("HMP_DeblurHigh.tsv", header=F, sep="\t")
Unoise_high_HMP <- read.table("HMP_UnoiseHigh.tsv", header=F, sep="\t")
HMP_High <- list(Dada_high_HMP, Deblur_high_HMP, Unoise_high_HMP)
HMP_High_plot <- plot_prop_curve(HMP_High, Pipeline_name_list, ylim=50)
HMP_High_plot

#Mock12

Dada_high_M12 <- read.table("Mock12_DadaHigh.tsv", header=F, sep="\t")
Deblur_high_M12 <- read.table("Mock12_DeblurHigh.tsv", header=F, sep="\t")
Unoise_high_M12 <- read.table("Mock12_UnoiseHigh.tsv", header=F, sep="\t")
M12_High <- list(Dada_high_M12, Deblur_high_M12, Unoise_high_M12)
M12_high_plot <- plot_prop_curve(M12_High, Pipeline_name_list, 50)
M12_high_plot

#not gonna bother with Mock9 due to the expected sequence list not being optimal?


#Zymock
Zy_Dada_High <- read.table("Zymock_DadaHigh.tsv", sep="\t")
Zy_Deblur_High <- read.table("Zymock_DeblurHigh.tsv", sep="\t")
Zy_Unoise_High <- read.table("Zymock_UnoiseHigh.tsv", sep="\t")
Zy_High <- list(Zy_Dada_High, Zy_Deblur_High, Zy_Unoise_High)
Zy_High_plot <- plot_prop_curve(Zy_High, Pipeline_name_list, 50)
Zy_High_plot



################################### Medium plots

#HMP MED

HMP_Dada_Med <- read.table("HMP_DadaMed.tsv", sep="\t")
HMP_Deblur_Med <- read.table("HMP_DeblurMed.tsv", sep="\t")
HMP_Unoise_Med <- read.table("HMP_UnoiseMed.tsv", sep="\t")
HMP_Open_Med <- read.table("HMP_OpenMed.tsv", sep="\t")
HMP_Med <- list(HMP_Dada_Med, HMP_Deblur_Med, HMP_Unoise_Med, HMP_Open_Med)
HMP_Med_Plot <- plot_prop_curve(HMP_Med, Pipeline_name_list, ylim=c(175,175,175,175))
HMP_Med_Plot


# Zymock Med
Zy_Med_plot_rec <- recordPlot(Zy_Med_plot)
Zy_Dada_Med <- read.table("Zymock_DadaMed.tsv", sep="\t")
Zy_Deblur_Med <- read.table("Zymock_DeblurMed.tsv", sep="\t")
Zy_Unoise_Med <- read.table("Zymock_UnoiseMed.tsv", sep="\t")
Zy_Open_Med <- read.table("Zymock_OpenMed.tsv", sep="\t")
Zy_Med <- list(Zy_Dada_Med, Zy_Deblur_Med, Zy_Unoise_High, Zy_Open_Med)
Zy_Med_plot <- plot_prop_curve(Zy_Med, Pipeline_name_list, ylim=c(75,75,75,75))
Zy_Med_plot


# Mock12 Med

M12_Dada_Med <- read.table("Mock12_DadaMed.tsv", sep="\t")
M12_Deblur_Med <- read.table("Mock12_DeblurMed.tsv", sep="\t")
M12_Unoise_Med <- read.table("Mock12_UnoiseMed.tsv", sep="\t")
M12_Open_Med <- read.table("Mock12_OpenMed.tsv", sep="\t")
M12_Med <- list(M12_Dada_Med, M12_Deblur_Med, M12_Unoise_Med, M12_Open_Med)
M12_Med_plot <- plot_prop_curve(M12_Med, Pipeline_name_list, ylim=c(45,45,45,4500))
M12_Med_plot

#Combined all Med plots together?


#seems really messy...
Master_Med_plot <- plot_grid(plotlist=list(M12_Med_plot, Zy_Med_plot, HMP_Med_Plot))
Master_Med_plot

################################### Low plots

#HMP Low
Low_HMP_Dada <- read.table("HMP_DadaLow.tsv", sep="\t")
Low_HMP_Deblur <- read.table("HMP_DeblurLow.tsv", sep="\t")
Low_HMP_Unoise <- read.table("HMP_UnoiseLow.tsv", sep="\t")
Low_HMP <- list(Low_HMP_Dada, Low_HMP_Deblur, Low_HMP_Unoise)
Low_HMP_plot <- plot_prop_curve(Low_HMP, Pipeline_name_list, 100)
Low_HMP_plot


#Zymock Low
Low_Zy_Dada <- read.table("Zymock_DadaLow.tsv", sep="\t")
Low_Zy_Deblur <- read.table("Zymock_DeblurLow.tsv", sep="\t")
Low_Zy_Unoise <- read.table("Zymock_UnoiseLow.tsv", sep="\t")
Low_Zy <- list(Low_Zy_Dada, Low_Zy_Deblur, Low_Zy_Unoise)
Low_Zy_plot <- plot_prop_curve(Low_Zy, Pipeline_name_list, 100)
Low_Zy_plot



####################################### Done