library(reshape2)
library(ggplot2)
library(grid)
library(cowplot)
library(gridExtra)
library(plyr)
files <- list.files(path = "~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/", pattern = "per_*")
setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")

HMP_DADA_HIGH <- makeFullTable("per_HMP_Dada_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_HIGH <- makeFullTable("per_HMP_Deblur_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_HIGH <- makeFullTable("per_HMP_Unoise_High_Org_col.tsv", "../expected/HMP_Expected.tsv")

StackedBarplotForHMP(HMP_DADA_HIGH, HMP_DEBLUR_HIGH, HMP_UNOISE_HIGH, "HMP_High")

HMP_DADA_MED <- makeFullTable("per_HMP_Dada_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_MED <- makeFullTable("per_HMP_Deblur_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_MED <- makeFullTable("per_HMP_Unoise_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")

StackedBarplotForHMP(HMP_DADA_MED, HMP_DEBLUR_MED, HMP_UNOISE_MED, "HMP_Med")

HMP_DADA_LOW <- makeFullTable("per_HMP_Dada_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_LOW <- makeFullTable("per_HMP_Deblur_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_LOW <- makeFullTable("per_HMP_Unoise_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")

StackedBarplotForHMP(HMP_DADA_LOW, HMP_DEBLUR_LOW, HMP_UNOISE_LOW, "HMP_Low")

Zymock_DADA_HIGH <- makeFullTable("per_Zymock_Dada_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_HIGH <- makeFullTable("per_Zymock_Deblur_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_HIGH <- makeFullTable("per_Zymock_Unoise_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")

StackedBarplotForZymock(Zymock_DADA_HIGH, Zymock_DEBLUR_HIGH, Zymock_UNOISE_HIGH, "Zymock_High")

Zymock_DADA_MED <- makeFullTable("per_Zymock_Dada_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_MED <- makeFullTable("per_Zymock_Deblur_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_MED <- makeFullTable("per_Zymock_Unoise_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")

StackedBarplotForZymock(Zymock_DADA_MED, Zymock_DEBLUR_MED, Zymock_UNOISE_MED, "Zymock_Med")

Zymock_DADA_LOW <- makeFullTable("per_Zymock_Dada_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_LOW <- makeFullTable("per_Zymock_Deblur_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_LOW <- makeFullTable("per_Zymock_Unoise_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")

StackedBarplotForZymock(Zymock_DADA_LOW, Zymock_DEBLUR_LOW, Zymock_UNOISE_LOW, "Zymock_Low")

Mock9_DADA_HIGH <- makeFullTable("per_Mock-9_Dada_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_HIGH <- makeFullTable("per_Mock-9_Deblur_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_HIGH <- makeFullTable("per_Mock-9_Unoise_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")

StackedBarplotForZymock(Mock9_DADA_HIGH, Mock9_DEBLUR_HIGH, Mock9_UNOISE_HIGH, "Mock-9_High")

Mock9_DADA_MED <- makeFullTable("per_Mock-9_Dada_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_MED <- makeFullTable("per_Mock-9_Deblur_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_MED <- makeFullTable("per_Mock-9_Unoise_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")

StackedBarplotForZymock(Mock9_DADA_MED, Mock9_DEBLUR_MED, Mock9_UNOISE_MED, "Mock-9_Med")

Mock9_DADA_LOW <- makeFullTable("per_Mock-9_Dada_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_LOW <- makeFullTable("per_Mock-9_Deblur_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_LOW <- makeFullTable("per_Mock-9_Unoise_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")

StackedBarplotForZymock(Mock9_DADA_LOW, Mock9_DEBLUR_LOW, Mock9_UNOISE_LOW, "Mock-9_Low")

Mock12_DADA_HIGH <- makeFullTable("per_Mock-12_Dada_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_HIGH <- makeFullTable("per_Mock-12_Deblur_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_HIGH <- makeFullTable("per_Mock-12_Unoise_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")

StackedBarplotForMock12(Mock12_DADA_HIGH, Mock12_DEBLUR_HIGH, Mock12_UNOSIE_HIGH, "Mock-12_High")

Mock12_DADA_MED <- makeFullTable("per_Mock-12_Dada_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_MED <- makeFullTable("per_Mock-12_Deblur_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_MED <- makeFullTable("per_Mock-12_Unoise_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")

StackedBarplotForMock12(Mock12_DADA_MED, Mock12_DEBLUR_MED, Mock12_UNOSIE_MED, "Mock-12_Med")

Mock12_DADA_LOW <- makeFullTable("per_Mock-12_Dada_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_LOW <- makeFullTable("per_Mock-12_Deblur_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_LOW <- makeFullTable("per_Mock-12_Unoise_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")

StackedBarplotForMock12(Mock12_DADA_LOW, Mock12_DEBLUR_LOW, Mock12_UNOSIE_LOW, "Mock-12_Low")

#combines the sample table with the expected table
#use this to do sqrt((expected-observed)^2)
makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  target <- expected$Organism
  final <- final[match(target, final$Organism),]
  rownames(final) <- NULL
  return(final)
}


StackedBarplotForHMP <- function(tab1,tab2,tab3, filt){
  qual_col <- c("#a6cee3",
                "#1f78b4",
                "#b2df8a",
                "#33a02c",
               "#fb9a99",
               "#e31a1c",
                "#fdbf6f",
                "#ff7f00",
                "#cab2d6",
               "#6a3d9a",
               "#ffff99",
               "#b15928")
               
  
  qual_col_repeat <- c(qual_col, qual_col[1:9])
  
  #table 1
  tab1_filt <- tab1[-grep("Total", tab1$Organism),]
  tab1_melt <- melt(data=tab1_filt, measure.vars=c("Sample1", "Sample2", "Sample3", "Sample4", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot1 <- ggplot(tab1_melt, aes(x = sample, 
                             y = rel_abun, 
                             fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
   scale_fill_manual(values=qual_col_repeat)
  
  #table 2
  tab2_filt <- tab2[-grep("Total", tab2$Organism),]
  tab2_melt <- melt(data=tab2_filt, measure.vars=c("Sample1", "Sample2", "Sample3", "Sample4", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot2 <- ggplot(tab2_melt, aes(x = sample, 
                            y = rel_abun, 
                             fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #table 3
  tab3_filt <- tab3[-grep("Total", tab3$Organism),]
  tab3_melt <- melt(data=tab3_filt, measure.vars=c("Sample1", "Sample2", "Sample3", "Sample4", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot3 <- ggplot(tab3_melt, aes(x = sample, 
                              y = rel_abun, 
                              fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #legend
  dumby <- ggplot(tab3_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=qual_col_repeat)
  
  legend <- get_legend(dumby)
  
  grid_plot <- plot_grid(plot1, plot2, plot3,legend, labels = c("A", "B", "C", filt))
  grid_plot
  save_plot(paste(filt, "png", sep ="."), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 1.3)
}


StackedBarplotForZymock <- function(tab1,tab2,tab3, filt)
  {
  qual_col <- c("#a6cee3",
                "#1f78b4",
                "#b2df8a",
                "#33a02c",
                "#fb9a99",
                "#e31a1c",
                "#fdbf6f",
                "#ff7f00",
                "#cab2d6",
                "#6a3d9a",
                "#ffff99",
                "#b15928")
  
  qual_col_repeat <- c(qual_col, qual_col[1:9])
  
  #table 1
  tab1_filt <- tab1[-grep("Total", tab1$Organism),]
  tab1_melt <- melt(data=tab1_filt, measure.vars=c("Sample1", "Sample2", "Sample3", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot1 <- ggplot(tab1_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #table 2
  tab2_filt <- tab2[-grep("Total", tab2$Organism),]
  tab2_melt <- melt(data=tab2_filt, measure.vars=c("Sample1", "Sample2", "Sample3", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot2 <- ggplot(tab2_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #table 3
  tab3_filt <- tab3[-grep("Total", tab3$Organism),]
  tab3_melt <- melt(data=tab3_filt, measure.vars=c("Sample1", "Sample2", "Sample3", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot3 <- ggplot(tab3_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #legend
  dumby <- ggplot(tab3_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=qual_col_repeat)
  
  legend <- get_legend(dumby)
  
  grid_plot <- plot_grid(plot1, plot2, plot3,legend, labels = c("A", "B", "C", filt))
  grid_plot
  save_plot(paste(filt, "png", sep ="."), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 1.3)
}

StackedBarplotForMock12 <- function(tab1,tab2,tab3, filt)
{
  qual_col <- c("#a6cee3",
                "#1f78b4",
                "#b2df8a",
                "#33a02c",
                "#fb9a99",
                "#e31a1c",
                "#fdbf6f",
                "#ff7f00",
                "#cab2d6",
                "#6a3d9a",
                "#ffff99",
                "#b15928")
  
  qual_col_repeat <- c(qual_col, qual_col[1:9], qual_col[1:9])
  
  #table 1
  tab1_filt <- tab1[-grep("Total", tab1$Organism),]
  tab1_melt <- melt(data=tab1_filt, measure.vars=c("Sample1", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot1 <- ggplot(tab1_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #table 2
  tab2_filt <- tab2[-grep("Total", tab2$Organism),]
  tab2_melt <- melt(data=tab2_filt, measure.vars=c("Sample1", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot2 <- ggplot(tab2_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #table 3
  tab3_filt <- tab3[-grep("Total", tab3$Organism),]
  tab3_melt <- melt(data=tab3_filt, measure.vars=c("Sample1", "Expected"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot3 <- ggplot(tab3_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity", show.legend = F) +
    scale_fill_manual(values=qual_col_repeat)
  
  #legend
  dumby <- ggplot(tab3_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=qual_col_repeat)
  
  legend <- get_legend(dumby)
  
  grid_plot <- plot_grid(plot1, plot2, plot3, legend, scale = c(1, 1, 1, 0.1), labels = c("A", "B", "C", filt))
  grid_plot
  save_plot(paste(filt, "png", sep ="."), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 2)
}






