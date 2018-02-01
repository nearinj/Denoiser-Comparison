#script that gets mean of each sample and graphs all them in one graph..
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)
library(plyr)

setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")

Mock9makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  final[is.na(final)] <- 0
  return(final)
}

makeFullTable <- function(add_t,add_E){
  table <- read.table(paste(add_t), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expected <- read.table(paste(add_E), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  final <- merge(table, expected, all = TRUE)
  final[is.na(final)] <- 0
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
    Dada$Organism <- gsub("_", " ", Dada$Organism)
    Dada$Organism <- factor(Dada$Organism, levels= Dada$Organism)
    colnames(Dada)[c(4, 6)] <- c("DADA2", "UNOISE3")
    
  return(Dada)
}

Mock9mergePipes <- function(Dada, Deblur, Unoise) {
  colnames(Dada)[ncol(Dada)] <- "DADA2"
  colnames(Deblur)[ncol(Deblur)] <- "Deblur"
  colnames(Unoise)[ncol(Unoise)] <- "UNOISE3"
  Dada$Otus <- NULL
  Deblur$Otus <- NULL
  Unoise$Otus <- NULL
  final <- rbind.fill(Dada, Deblur, Unoise)
  final[is.na(final)] <- 0
  final <- ddply(final, .(Organism), summarize, DADA2=sum(DADA2), Deblur=sum(Deblur), UNOISE3=sum(UNOISE3), Expected=sum(Expected)/3)
  final[3,1] <- "Galactomyces geotrichum Type 1"
  final[3, "Expected"] <- 0.0625
  final[4,1] <- "Galactomyces geotrichum Type 2"
  final[4, "Expected"] <- 0.0625
  final[5,1] <- "Galactomyces geotrichum Type 3"
  final[5, "Expected"] <- 0.0625
  final[6,1] <- "Galactomyces geotrichum Type 4"
  final[6, "Expected"] <- 0.0625
  final[12,1] <- "Hyphopichia_burtonii"
  final[12, "Expected"] <- 0.0625
  final <- final[-c(10, 11, 8),]
  final$Organism <- gsub("_", " ", final$Organism)
  target <- final[final$Organism != "Non Reference", 1]
  final$Organism <- factor(final$Organism, levels=c("Non Reference", target))
  return(final)
}

collapseStaph <- function(tab) {
  New_Staph <- tab[17, c("Expected", "DADA2", "Deblur", "UNOISE3")] +tab[18, c("Expected", "DADA2", "Deblur", "UNOISE3") ]
  tab <- tab[-c(17,18), ]
  New_Staph$Organism <- "S.aureus/epidermidis"
  New_Staph$Otus <- "N/a"
  tab <- merge(tab, New_Staph, all=T)
}


HMP_DADA_HIGH <- makeFullTable("per_HMP_Dada_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_HIGH_MEAN <- getMeanHMP(HMP_DADA_HIGH)
HMP_DEBLUR_HIGH <- makeFullTable("per_HMP_Deblur_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_HIGH_MEAN <- getMeanHMP(HMP_DADA_HIGH)
HMP_UNOISE_HIGH <- makeFullTable("per_HMP_Unoise_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_HIGH_MEAN <- getMeanHMP(HMP_UNOISE_HIGH)

HMP_HIGH <- mergePipes(HMP_DADA_HIGH_MEAN, HMP_DEBLUR_HIGH_MEAN, HMP_UNOISE_HIGH_MEAN)
HMP_HIGH <- collapseStaph(HMP_HIGH)

HMP_DADA_MED <- makeFullTable("per_HMP_Dada_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_MED_MEAN <- getMeanHMP(HMP_DADA_MED)
HMP_DEBLUR_MED <- makeFullTable("per_HMP_Deblur_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_MED_MEAN <- getMeanHMP(HMP_DEBLUR_MED)
HMP_UNOISE_MED <- makeFullTable("per_HMP_Unoise_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_MED_MEAN <- getMeanHMP(HMP_UNOISE_MED)

HMP_MED <- mergePipes(HMP_DADA_MED_MEAN, HMP_DEBLUR_MED_MEAN, HMP_UNOISE_MED_MEAN)
HMP_MED <- collapseStaph(HMP_MED)

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

Mock9_DADA_HIGH <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Dada_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_HIGH_MEAN <- getMeanZymock(Mock9_DADA_HIGH)
Mock9_DEBLUR_HIGH <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Deblur_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_HIGH_MEAN <- getMeanZymock(Mock9_DEBLUR_HIGH)
Mock9_UNOISE_HIGH <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Unoise_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_HIGH_MEAN <- getMeanZymock(Mock9_UNOISE_HIGH)

Mock9_HIGH <- Mock9mergePipes(Mock9_DADA_HIGH_MEAN, Mock9_DEBLUR_HIGH_MEAN, Mock9_UNOISE_HIGH_MEAN)

Mock12_DADA_HIGH <- makeFullTable("per_Mock-12_Dada_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_HIGH_MEAN <- getMeanMock12(Mock12_DADA_HIGH)
Mock12_DEBLUR_HIGH <- makeFullTable("per_Mock-12_Deblur_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_HIGH_MEAN <- getMeanMock12(Mock12_DEBLUR_HIGH)
Mock12_UNOSIE_HIGH <- makeFullTable("per_Mock-12_Unoise_High_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_HIGH_MEAN <- getMeanMock12(Mock12_UNOSIE_HIGH)

Mock12_HIGH <- mergePipes(Mock12_DADA_HIGH_MEAN, Mock12_DEBLUR_HIGH_MEAN, Mock12_UNOSIE_HIGH_MEAN)

Zymock_DADA_MED <- makeFullTable("per_Zymock_Dada_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
ZYMOCK_DADA_MED_MEAN <- getMeanZymock(Zymock_DADA_MED)
Zymock_DEBLUR_MED <- makeFullTable("per_Zymock_Deblur_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_MED_MEAN <- getMeanZymock(Zymock_DEBLUR_MED)
Zymock_UNOISE_MED <- makeFullTable("per_Zymock_Unoise_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_MED_MEAN <- getMeanZymock(Zymock_UNOISE_MED)

Zymock_MED <- mergePipes(ZYMOCK_DADA_MED_MEAN, Zymock_DEBLUR_MED_MEAN, Zymock_UNOISE_MED_MEAN)

Mock9_DADA_MED <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Dada_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_MED_MEAN <- getMeanZymock(Mock9_DADA_MED)
Mock9_DEBLUR_MED <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Deblur_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_MED_MEAN <- getMeanZymock(Mock9_DEBLUR_MED)
Mock9_UNOISE_MED <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Unoise_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_MED_MEAN <- getMeanZymock(Mock9_UNOISE_MED)

Mock9_MED <- Mock9mergePipes(Mock9_DADA_MED_MEAN, Mock9_DEBLUR_MED_MEAN, Mock9_UNOISE_MED_MEAN)

Mock12_DADA_MED <- makeFullTable("per_Mock-12_Dada_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_MED_MEAN <- getMeanMock12(Mock12_DADA_MED)
Mock12_DEBLUR_MED <- makeFullTable("per_Mock-12_Deblur_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_MED_MEAN <- getMeanMock12(Mock12_DEBLUR_MED)
Mock12_UNOSIE_MED <- makeFullTable("per_Mock-12_Unoise_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_MED_MEAN <- getMeanMock12(Mock12_UNOSIE_MED)

Mock12_MED <- mergePipes(Mock12_DADA_MED_MEAN, Mock12_DEBLUR_MED_MEAN, Mock12_UNOSIE_MED_MEAN)

Zymock_DADA_LOW <- makeFullTable("per_Zymock_Dada_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
ZYMOCK_DADA_LOW_MEAN <- getMeanZymock(Zymock_DADA_LOW)
Zymock_DEBLUR_LOW <- makeFullTable("per_Zymock_Deblur_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_LOW_MEAN <- getMeanZymock(Zymock_DEBLUR_LOW)
Zymock_UNOISE_LOW <- makeFullTable("per_Zymock_Unoise_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_LOW_MEAN <- getMeanZymock(Zymock_UNOISE_LOW)

Zymock_LOW <- mergePipes(ZYMOCK_DADA_LOW_MEAN, Zymock_DEBLUR_LOW_MEAN, Zymock_UNOISE_LOW_MEAN)

Mock9_DADA_LOW <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Dada_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_LOW_MEAN <- getMeanZymock(Mock9_DADA_LOW)
Mock9_DEBLUR_LOW <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Deblur_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_LOW_MEAN <- getMeanZymock(Mock9_DEBLUR_LOW)
Mock9_UNOISE_LOW <- Mock9makeFullTable("../incl-db-collapsed/Mock-9_Unoise_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_LOW_MEAN <- getMeanZymock(Mock9_UNOISE_LOW)

Mock9_LOW <- Mock9mergePipes(Mock9_DADA_LOW_MEAN, Mock9_DEBLUR_LOW_MEAN, Mock9_UNOISE_LOW_MEAN)

Mock12_DADA_LOW <- makeFullTable("per_Mock-12_Dada_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_LOW_MEAN <- getMeanMock12(Mock12_DADA_LOW)
Mock12_DEBLUR_LOW <- makeFullTable("per_Mock-12_Deblur_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_LOW_MEAN <- getMeanMock12(Mock12_DEBLUR_LOW)
Mock12_UNOSIE_LOW <- makeFullTable("per_Mock-12_Unoise_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_LOW_MEAN <- getMeanMock12(Mock12_UNOSIE_LOW)

Mock12_LOW <- mergePipes(Mock12_DADA_LOW_MEAN, Mock12_DEBLUR_LOW_MEAN, Mock12_UNOSIE_LOW_MEAN)

Medium <- StackedBarplot(HMP_MED, Zymock_MED, Mock9_MED, Mock12_MED, "Med")
Medium

High <- StackedBarplot(HMP_HIGH, Zymock_HIGH, Mock9_HIGH, Mock12_HIGH, "High")
High

Low <- StackedBarplot(HMP_LOW, Zymock_LOW, Mock9_LOW, Mock12_LOW, "Low")
Low


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
  HMPTable_melt <- melt(data=HMPTable_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot1 <- ggplot(HMPTable_melt, aes(x = sample, 
                                 y = rel_abun, 
                                 fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=12)) +
                                                      labs(x="Pipelines", y="Relative Abundance")
  
  plot1
  
  
  #table 2 Zymock_Table
  ZymockTable_filt <- ZymockTable[-grep("Total", ZymockTable$Organism),]
  ZymockTable_melt <- melt(ZymockTable_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3"), 
                        id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot2 <- ggplot(ZymockTable_melt, aes(x = sample, 
                                     y = rel_abun, 
                                     fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=12)) + 
                                                      labs(x="Pipelines", y="Relative Abundance")
  
  plot2
  
  #table 3 Mock9
  Mock9Table_filt <- Mock9Table[-grep("Total", Mock9Table$Organism),]
  
  Mock9Table_melt <- melt(data=Mock9Table_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3"), 
                          id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot3 <- ggplot(Mock9Table_melt, aes(x = sample, 
                                       y = rel_abun, 
                                       fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=12)) +
                                                      labs(x = "Pipelines", y="Relative Abundance")
  
  plot3
  
  #table 4 Mock12
  Mock12Table_filt <- Mock12Table[-grep("Total", Mock12Table$Organism),]
  Mock12Table_melt <- melt(data=Mock12Table_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3"), 
                           id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot4 <- ggplot(Mock12Table_melt, aes(x = sample, 
                                        y = rel_abun, 
                                        fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                      legend.title = element_text(size=10), axis.text.x = element_text(size=12)) + 
                                                      labs(x = "Pipelines", y="Relative Abundance")
  
  plot4
  
  
  
  grid_plot <- plot_grid(plot1, plot4, plot3, plot2, labels=c("A","B","C","D"))
  grid_plot
  #save_plot(paste("New",filt, ".png", sep =""), grid_plot, ncol = 2, nrow = 2, base_aspect_ratio = 2)
  return(grid_plot)
}



