#script that gets mean of each sample and graphs all them in one graph..
library(ggplot2)
library(scales)
library(reshape2)
library(cowplot)
library(plyr)

setwd("~/projects/DenoiseCompare/Data_Analysis/fixed_ASV_analysis/Abundance/collapsed/")

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

mergePipes_Med <- function(Dada, Deblur, Unoise, Open) {
  colnames(Dada)[ncol(Dada)] <- "Dada"
  Dada$Deblur <- Deblur$mean
  Dada$Unoise <- Unoise$mean
  Dada$Open <- Open$mean
  Dada <- replace(Dada, is.na(Dada), 0)
  Dada$Organism <- gsub("_", " ", Dada$Organism)
  Dada$Organism <- factor(Dada$Organism, levels= Dada$Organism)
  colnames(Dada)[c(4, 6)] <- c("DADA2", "UNOISE3")
  
  return(Dada)
}

Mock9mergePipes_Med <- function(Dada, Deblur, Unoise, Open) {
  colnames(Dada)[ncol(Dada)] <- "DADA2"
  colnames(Deblur)[ncol(Deblur)] <- "Deblur"
  colnames(Unoise)[ncol(Unoise)] <- "UNOISE3"
  colnames(Open)[ncol(Open)] <- "Open"
  Dada$Otus <- NULL
  Deblur$Otus <- NULL
  Unoise$Otus <- NULL
  final <- rbind.fill(Dada, Deblur, Unoise, Open)
  final[is.na(final)] <- 0
  final <- ddply(final, .(Organism), summarize, DADA2=sum(DADA2), Deblur=sum(Deblur), UNOISE3=sum(UNOISE3), 
                 Open=sum(Open), Expected=sum(Expected)/4)
  
  #remove anything that matches Fungi sp|K
  
  
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
  #remove old expected and Emeercella group with 5 mapping reads
  final <- final[-c(10, 11, 8),]
  final$Organism <- gsub("_", " ", final$Organism)
  target <- final[final$Organism != "Non Reference", 1]
  final$Organism <- factor(final$Organism, levels=c("Non Reference", target))
  return(final)
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
  #remove old expected and Emeercella group with 5 mapping reads
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


collapseStaph_Med <- function(tab) {
  New_Staph <- tab[17, c("Expected", "DADA2", "Deblur", "UNOISE3", "Open")] +tab[18, c("Expected", "DADA2", "Deblur", "UNOISE3", "Open") ]
  tab <- tab[-c(17,18), ]
  New_Staph$Organism <- "S.aureus/epidermidis"
  New_Staph$Otus <- "N/a"
  tab <- merge(tab, New_Staph, all=T)
}

StackedBarplot <- function(HMPTable, ZymockTable, Mock9Table, Mock12Table, filt ) {
  qual_col <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
  
  
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

StackedBarplot_MED <- function(HMPTable, ZymockTable, Mock9Table, Mock12Table, filt ) {
  qual_col <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
  
  
  qual_col_repeat <- c(qual_col, qual_col, qual_col)
  
  #table 1 HMP Table
  HMPTable_filt <- HMPTable[-grep("Total", HMPTable$Organism),]
  HMPTable_melt <- melt(data=HMPTable_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open"), 
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
  ZymockTable_melt <- melt(ZymockTable_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open"), 
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
  
  Mock9Table_melt <- melt(data=Mock9Table_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open"), 
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
  Mock12Table_melt <- melt(data=Mock12Table_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open"), 
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

Mock9_DADA_HIGH <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Dada_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_HIGH_MEAN <- getMeanZymock(Mock9_DADA_HIGH)
Mock9_DEBLUR_HIGH <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Deblur_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_HIGH_MEAN <- getMeanZymock(Mock9_DEBLUR_HIGH)
Mock9_UNOISE_HIGH <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Unoise_High_Org_col.tsv", "../expected/Mock9_Expected.tsv")
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


Mock9_DADA_MED <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Dada_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_MED_MEAN <- getMeanZymock(Mock9_DADA_MED)
Mock9_DEBLUR_MED <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Deblur_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_MED_MEAN <- getMeanZymock(Mock9_DEBLUR_MED)
Mock9_UNOISE_MED <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Unoise_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_MED_MEAN <- getMeanZymock(Mock9_UNOISE_MED)



Mock12_DADA_MED <- makeFullTable("per_Mock-12_Dada_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_MED_MEAN <- getMeanMock12(Mock12_DADA_MED)
Mock12_DEBLUR_MED <- makeFullTable("per_Mock-12_Deblur_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_MED_MEAN <- getMeanMock12(Mock12_DEBLUR_MED)
Mock12_UNOSIE_MED <- makeFullTable("per_Mock-12_Unoise_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_MED_MEAN <- getMeanMock12(Mock12_UNOSIE_MED)


Zymock_DADA_LOW <- makeFullTable("per_Zymock_Dada_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
ZYMOCK_DADA_LOW_MEAN <- getMeanZymock(Zymock_DADA_LOW)
Zymock_DEBLUR_LOW <- makeFullTable("per_Zymock_Deblur_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_DEBLUR_LOW_MEAN <- getMeanZymock(Zymock_DEBLUR_LOW)
Zymock_UNOISE_LOW <- makeFullTable("per_Zymock_Unoise_Low_Org_col.tsv", "../expected/Zymock_Expected.tsv")
Zymock_UNOISE_LOW_MEAN <- getMeanZymock(Zymock_UNOISE_LOW)

Zymock_LOW <- mergePipes(ZYMOCK_DADA_LOW_MEAN, Zymock_DEBLUR_LOW_MEAN, Zymock_UNOISE_LOW_MEAN)

Mock9_DADA_LOW <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Dada_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DADA_LOW_MEAN <- getMeanZymock(Mock9_DADA_LOW)
Mock9_DEBLUR_LOW <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Deblur_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_DEBLUR_LOW_MEAN <- getMeanZymock(Mock9_DEBLUR_LOW)
Mock9_UNOISE_LOW <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Unoise_Low_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_UNOISE_LOW_MEAN <- getMeanZymock(Mock9_UNOISE_LOW)

Mock9_LOW <- Mock9mergePipes(Mock9_DADA_LOW_MEAN, Mock9_DEBLUR_LOW_MEAN, Mock9_UNOISE_LOW_MEAN)

Mock12_DADA_LOW <- makeFullTable("per_Mock-12_Dada_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DADA_LOW_MEAN <- getMeanMock12(Mock12_DADA_LOW)
Mock12_DEBLUR_LOW <- makeFullTable("per_Mock-12_Deblur_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_DEBLUR_LOW_MEAN <- getMeanMock12(Mock12_DEBLUR_LOW)
Mock12_UNOSIE_LOW <- makeFullTable("per_Mock-12_Unoise_Low_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_UNOSIE_LOW_MEAN <- getMeanMock12(Mock12_UNOSIE_LOW)

Mock12_LOW <- mergePipes(Mock12_DADA_LOW_MEAN, Mock12_DEBLUR_LOW_MEAN, Mock12_UNOSIE_LOW_MEAN)

#Open Tables
HMP_OPEN_MED <- makeFullTable("per_HMP_Open_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_OPEN_MED_Mean <- getMeanHMP(HMP_OPEN_MED)

Mock12_OPEN_MED <- makeFullTable("per_Mock-12_Open_Med_Org_col.tsv", "../expected/Mock-12_Expected.tsv")
Mock12_OPEN_MED_Mean <- getMeanMock12(Mock12_OPEN_MED)

Mock9_OPEN_MED <- Mock9makeFullTable("incl-db-collapsed/Mock-9_Open_Med_Org_col.tsv", "../expected/Mock9_Expected.tsv")
Mock9_OPEN_MED_MEAN <- getMeanZymock(Mock9_OPEN_MED)

ZYMOCK_OPEN_MED <- makeFullTable("per_Zymock_Open_Med_Org_col.tsv", "../expected/Zymock_Expected.tsv")
ZYMOCK_OPEN_MED_MEAN <- getMeanZymock(ZYMOCK_OPEN_MED)


#Merge the Medium filter pipelines that include open-ref 
Mock9_MED <- Mock9mergePipes_Med(Mock9_DADA_MED_MEAN, Mock9_DEBLUR_MED_MEAN, Mock9_UNOISE_MED_MEAN, Mock9_OPEN_MED_MEAN)
Mock12_MED <- mergePipes_Med(Mock12_DADA_MED_MEAN, Mock12_DEBLUR_MED_MEAN, Mock12_UNOSIE_MED_MEAN, Mock12_OPEN_MED_Mean)
Zymock_MED <- mergePipes_Med(ZYMOCK_DADA_MED_MEAN, Zymock_DEBLUR_MED_MEAN, Zymock_UNOISE_MED_MEAN, ZYMOCK_OPEN_MED_MEAN)
HMP_MED <- mergePipes_Med(HMP_DADA_MED_MEAN, HMP_DEBLUR_MED_MEAN, HMP_UNOISE_MED_MEAN, HMP_OPEN_MED_Mean)
HMP_MED <- collapseStaph_Med(HMP_MED)


Medium <- StackedBarplot_MED(HMP_MED, Zymock_MED, Mock9_MED, Mock12_MED, "Med")
Medium

High <- StackedBarplot(HMP_HIGH, Zymock_HIGH, Mock9_HIGH, Mock12_HIGH, "High")
High

Low <- StackedBarplot(HMP_LOW, Zymock_LOW, Mock9_LOW, Mock12_LOW, "Low")
Low


StackedBarplot <- function(HMPTable, ZymockTable, Mock9Table, Mock12Table, filt ) {
  qual_col <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                 "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
  
  
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

StackedBarplot_MED <- function(HMPTable, ZymockTable, Mock9Table, Mock12Table, filt ) {
  qual_col <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
  
  
  qual_col_repeat <- c(qual_col, qual_col, qual_col)
  
  #table 1 HMP Table
  HMPTable_filt <- HMPTable[-grep("Total", HMPTable$Organism),]
  HMPTable_melt <- melt(data=HMPTable_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open_Ref_OTU"), 
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
  ZymockTable_melt <- melt(ZymockTable_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open_Ref_OTU"), 
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
  
  Mock9Table_melt <- melt(data=Mock9Table_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open_Ref_OTU"), 
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
  Mock12Table_melt <- melt(data=Mock12Table_filt, measure.vars=c("Expected", "DADA2", "Deblur", "UNOISE3", "Open_Ref_OTU"), 
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






