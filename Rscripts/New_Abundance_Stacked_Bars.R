#script that gets mean of each sample and graphs all them in one graph..
library(ggplot2)

makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  target <- expected$Organism
  final <- final[match(target, final$Organism),]
  rownames(final) <- NULL
  return(final)
}

getMeanHMP <- function(table) {
  table$mean <- rowMeans(table[, c(3,4,5,6)])
  table <- table[, -c(3,4,5,6)]
  return(table)
}

getMeanZymock <- function(table) {
  table$mean <- rowMeans(table[, c(3,4,5)])
  table <- table[, -c(3,4,5)]
  return(table)
}

getMeanMock12 <- function(table) {
  table$mean <- rowMeans(table[, c(3,3)])
  table <- table[, -3]
  return(table)
}

mergePipes <- function(Dada, Deblur, Unoise) {
    colnames(Dada)[ncol(Dada)] <- "Dada"
    Dada$Deblur <- Deblur$mean
    Dada$Unoise <- Unoise$mean
    Dada <- replace(Dada, is.na(Dada), 0)
    Dada$Organism <- factor(Dada$Organism, levels= Dada$Organism)
  return(Dada)
}

setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")

HMP_DADA_HIGH <- makeFullTable("per_HMP_Dada_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_HIGH_MEAN <- getMeanHMP(HMP_DADA_HIGH)
HMP_DEBLUR_HIGH <- makeFullTable("per_HMP_Deblur_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_HIGH_MEAN <- getMeanHMP(HMP_DADA_HIGH)
HMP_UNOISE_HIGH <- makeFullTable("per_HMP_Unoise_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_HIGH_MEAN <- getMeanHMP(HMP_UNOISE_HIGH)

HMP_HIGH <- mergePipes(HMP_DADA_HIGH_MEAN, HMP_DEBLUR_HIGH_MEAN, HMP_UNOISE_HIGH_MEAN)

HMP_DADA_MED <- makeFullTable("per_HMP_Dada_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_MED_MEAN <- getMeanHMP(HMP_DADA_MED)
HMP_DEBLUR_MED <- makeFullTable("per_HMP_Deblur_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_MED_MEAN <- getMeanHMP(HMP_DEBLUR_MED)
HMP_UNOISE_MED <- makeFullTable("per_HMP_Unoise_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_MED_MEAN <- getMeanHMP(HMP_UNOISE_MED)

HMP_MED <- mergePipes(HMP_DADA_MED_MEAN, HMP_DEBLUR_MED_MEAN, HMP_UNOISE_MED_MEAN)

HMP_DADA_LOW <- makeFullTable("per_HMP_Dada_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_LOW_MEAN <- getMeanHMP(HMP_DADA_LOW)
HMP_DEBLUR_LOW <- makeFullTable("per_HMP_Deblur_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_LOW_MEAN <- getMeanHMP(HMP_DADA_LOW)
HMP_UNOISE_LOW <- makeFullTable("per_HMP_Unoise_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_LOW_MEAN <- getMeanHMP(HMP_UNOISE_LOW)

HMP_LOW <- mergePipes(HMP_DADA_LOW_MEAN, HMP_DEBLUR_LOW_MEAN, HMP_UNOISE_LOW_MEAN)


Zymock_DADA_HIGH <- makeFullTable("per_Zymock_Dada_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
ZYMOCK_DADA_HIGH_MEAN <- getMeanZymock(Zymock_DADA_HIGH)
Zymock_DEBLUR_HIGH <- makeFullTable("per_Zymock_Deblur_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_HIGH_MEAN <- getMeanZymock(Zymock_DEBLUR_HIGH)
Zymock_UNOISE_HIGH <- makeFullTable("per_Zymock_Unoise_High_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_HIGH_MEAN <- getMeanZymock(Zymock_UNOISE_HIGH)

Zymock_HIGH <- mergePipes(ZYMOCK_DADA_HIGH_MEAN, Zymock_DEBLUR_HIGH_MEAN, Zymock_UNOISE_HIGH_MEAN)

Mock9_DADA_HIGH <- makeFullTable("per_Mock-9_Dada_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_HIGH_MEAN <- getMeanZymock(Mock9_DADA_HIGH)
Mock9_DEBLUR_HIGH <- makeFullTable("per_Mock-9_Deblur_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_HIGH_MEAN <- getMeanZymock(Mock9_DEBLUR_HIGH)
Mock9_UNOISE_HIGH <- makeFullTable("per_Mock-9_Unoise_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_HIGH_MEAN <- getMeanZymock(Mock9_UNOISE_HIGH)

Mock9_HIGH <- mergePipes(Mock9_DADA_HIGH_MEAN, Mock9_DEBLUR_HIGH_MEAN, Mock9_UNOISE_HIGH_MEAN)

Mock12_DADA_HIGH <- makeFullTable("per_Mock-12_Dada_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_HIGH_MEAN <- getMeanMock12(Mock12_DADA_HIGH)
Mock12_DEBLUR_HIGH <- makeFullTable("per_Mock-12_Deblur_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_HIGH_MEAN <- getMeanMock12(Mock12_DEBLUR_HIGH)
Mock12_UNOSIE_HIGH <- makeFullTable("per_Mock-12_Unoise_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_HIGH_MEAN <- getMeanMock12(Mock12_UNOSIE_HIGH)

Mock12_HIGH <- mergePipes(Mock12_DADA_HIGH_MEAN, Mock12_DEBLUR_HIGH_MEAN, Mock12_UNOSIE_HIGH_MEAN)

StackedBarplot(HMP_HIGH, Zymock_HIGH, Mock9_HIGH, Mock12_HIGH, "High")


Zymock_DADA_MED <- makeFullTable("per_Zymock_Dada_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
ZYMOCK_DADA_MED_MEAN <- getMeanZymock(Zymock_DADA_MED)
Zymock_DEBLUR_MED <- makeFullTable("per_Zymock_Deblur_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_MED_MEAN <- getMeanZymock(Zymock_DEBLUR_MED)
Zymock_UNOISE_MED <- makeFullTable("per_Zymock_Unoise_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_MED_MEAN <- getMeanZymock(Zymock_UNOISE_MED)

Zymock_MED <- mergePipes(ZYMOCK_DADA_MED_MEAN, Zymock_DEBLUR_MED_MEAN, Zymock_UNOISE_MED_MEAN)

Mock9_DADA_MED <- makeFullTable("per_Mock-9_Dada_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_MED_MEAN <- getMeanZymock(Mock9_DADA_MED)
Mock9_DEBLUR_MED <- makeFullTable("per_Mock-9_Deblur_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_MED_MEAN <- getMeanZymock(Mock9_DEBLUR_MED)
Mock9_UNOISE_MED <- makeFullTable("per_Mock-9_Unoise_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_MED_MEAN <- getMeanZymock(Mock9_UNOISE_MED)

Mock9_MED <- mergePipes(Mock9_DADA_MED_MEAN, Mock9_DEBLUR_MED_MEAN, Mock9_UNOISE_MED_MEAN)

Mock12_DADA_MED <- makeFullTable("per_Mock-12_Dada_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_MED_MEAN <- getMeanMock12(Mock12_DADA_MED)
Mock12_DEBLUR_MED <- makeFullTable("per_Mock-12_Deblur_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_MED_MEAN <- getMeanMock12(Mock12_DEBLUR_MED)
Mock12_UNOSIE_MED <- makeFullTable("per_Mock-12_Unoise_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_MED_MEAN <- getMeanMock12(Mock12_UNOSIE_MED)

Mock12_MED <- mergePipes(Mock12_DADA_MED_MEAN, Mock12_DEBLUR_MED_MEAN, Mock12_UNOSIE_MED_MEAN)

StackedBarplot(HMP_MED, Zymock_MED, Mock9_MED, Mock12_MED, "Med")

StackedBarplot <- function(HMPTable, ZymockTable, Mock9Table, Mock12Table, filt ) {
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
  
  
  qual_col_repeat <- c(qual_col, qual_col, qual_col)
  
  #table 1 HMP Table
  HMPTable_filt <- HMPTable[-grep("Total", HMPTable$Organism),]
  HMPTable_melt <- melt(data=HMPTable_filt, measure.vars=c("Expected", "Dada", "Deblur", "Unoise"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot1 <- ggplot(HMPTable_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#000000", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=6)) +
                                                      labs(x="Pipelines", y="Relative Abundance")
  
  plot1
  
  
  #table 2 Zymock_Table
  ZymockTable_filt <- ZymockTable[-grep("Total", ZymockTable$Organism),]
  ZymockTable_melt <- melt(data=ZymockTable_filt, measure.vars=c("Expected", "Dada", "Deblur", "Unoise"), 
                        id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot2 <- ggplot(ZymockTable_melt, aes(x = sample, 
                                     y = rel_abun, 
                                     fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#000000", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=6)) + 
                                                      labs(x="Pipelines", y="Relative Abundance")
  
  plot2
  
  #table 3 Mock9
  Mock9Table_filt <- Mock9Table[-grep("Total", Mock9Table$Organism),]
  Mock9Table_melt <- melt(data=Mock9Table_filt, measure.vars=c("Expected", "Dada", "Deblur", "Unoise"), 
                          id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot3 <- ggplot(Mock9Table_melt, aes(x = sample, 
                                       y = rel_abun, 
                                       fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#000000", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=6)) +
                                                      labs(x = "Pipelines", y="Relative Abundance")
  
  plot3
  
  #table 4 Mock12
  Mock12Table_filt <- Mock12Table[-grep("Total", Mock12Table$Organism),]
  Mock12Table_melt <- melt(data=Mock12Table_filt, measure.vars=c("Expected", "Dada", "Deblur", "Unoise"), 
                           id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot4 <- ggplot(Mock12Table_melt, aes(x = sample, 
                                        y = rel_abun, 
                                        fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#000000", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=6)) + 
                                                      labs(x = "Pipelines", y="Relative Abundance")
  
  plot4
  
  
  
  
  grid_plot <- plot_grid(plot1, plot2, plot3, plot4, labels=c("A","B","C","D"))
  grid_plot
  save_plot(paste("New",filt, ".png", sep =""), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 2)
}



