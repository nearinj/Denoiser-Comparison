library(ggplot2)
library(scales)

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
  table <- table[, -c(2,3)]
  return(table)
}

collapseStaph <- function(tab) {
  New_Staph <- tab[17, c("Expected", "Dada", "Deblur", "Unoise")] +HMP_HIGH[18, c("Expected", "Dada", "Deblur", "Unoise") ]
  tab <- tab[-c(17,18), ]
  New_Staph$Organism <- "S.aureus/epidermidis"
  New_Staph$Otus <- "N/a"
  tab <- merge(tab, New_Staph, all=T)
}

setwd("~/projects/DenoiseCompare/Data_Analysis/Abundance/collapsed/")

#Make all high tables
HMP_DADA_HIGH <- makeFullTable("per_HMP_Dada_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_HIGH_MEAN <- getMeanHMP(HMP_DADA_HIGH)
HMP_DEBLUR_HIGH <- makeFullTable("per_HMP_Deblur_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_HIGH_MEAN <- getMeanHMP(HMP_DADA_HIGH)
HMP_UNOISE_HIGH <- makeFullTable("per_HMP_Unoise_High_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_HIGH_MEAN <- getMeanHMP(HMP_UNOISE_HIGH)

#Make all Med tables
HMP_DADA_MED <- makeFullTable("per_HMP_Dada_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_MED_MEAN <- getMeanHMP(HMP_DADA_MED)
HMP_DEBLUR_MED <- makeFullTable("per_HMP_Deblur_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_MED_MEAN <- getMeanHMP(HMP_DEBLUR_MED)
HMP_UNOISE_MED <- makeFullTable("per_HMP_Unoise_Med_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_MED_MEAN <- getMeanHMP(HMP_UNOISE_MED)

#Make all Low tables
HMP_DADA_LOW <- makeFullTable("per_HMP_Dada_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DADA_LOW_MEAN <- getMeanHMP(HMP_DADA_LOW)
HMP_DEBLUR_LOW <- makeFullTable("per_HMP_Deblur_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_DEBLUR_LOW_MEAN <- getMeanHMP(HMP_DEBLUR_LOW)
HMP_UNOISE_LOW <- makeFullTable("per_HMP_Unoise_Low_Org_col.tsv", "../expected/HMP_Expected.tsv")
HMP_UNOISE_LOW_MEAN <- getMeanHMP(HMP_UNOISE_LOW)

HMP_DADA2 <- combineByPipeline(HMP_DADA_LOW_MEAN, HMP_DADA_MED_MEAN, HMP_DADA_HIGH_MEAN) 
HMP_Deblur <- combineByPipeline(HMP_DEBLUR_LOW_MEAN, HMP_DEBLUR_MED_MEAN, HMP_DEBLUR_HIGH_MEAN)
HMP_UNOISE <- combineByPipeline(HMP_UNOISE_LOW_MEAN, HMP_UNOISE_MED_MEAN, HMP_UNOISE_HIGH_MEAN)

PlotByPipeline(HMP_DADA2, HMP_Deblur, HMP_UNOISE)

combineByPipeline <- function(Lowtab, Medtab, Hightab){
  
  colnames(Hightab)[2] <- "High"
  Hightab$Med <- Medtab$mean
  Hightab$Low <- Lowtab$mean
  Hightab <- replace(Hightab, is.na(Hightab), 0)
  Hightab$Organism <- factor(Hightab$Organism, levels=Hightab$Organism)
  return(Hightab)
}

PlotByPipeline <- function(dada, deblur, unoise){
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
  dada_filt <- dada[-grep("Total", dada$Organism),]
  dada_melt <- melt(data=dada_filt, measure.vars=c("High", "Med", "Low"), 
                        id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot1 <- ggplot(dada_melt, aes(x = sample, 
                                     y = rel_abun, 
                                     fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                                    legend.title = element_text(size=10), axis.text.x = element_text(size=6)) +
    labs(x="Filter Stringency", y="Relative Abundance")
  
  plot1
  
  deblur_filt <- deblur[-grep("Total", deblur$Organism),]
  deblur_melt <- melt(data=deblur_filt, measure.vars=c("High", "Med", "Low"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot2 <- ggplot(deblur_melt, aes(x = sample, 
                                     y = rel_abun, 
                                     fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                                    legend.title = element_text(size=10), axis.text.x = element_text(size=6)) +
    labs(x="Filter Stringency", y="Relative Abundance")
  
  plot2
  
  unoise_filt <- unoise[-grep("Total", unoise$Organism),]
  unoise_melt <- melt(data=unoise_filt, measure.vars=c("High", "Med", "Low"), 
                    id.vars = c("Organism"), value.name="rel_abun", variable.name="sample")
  plot3 <- ggplot(unoise_melt, aes(x = sample, 
                                     y = rel_abun, 
                                     fill = Organism)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("#bdbdbd", qual_col_repeat)) + theme(legend.key.size = unit(.5, "line"), legend.text=element_text(size=6), 
                                                                    legend.title = element_text(size=10), axis.text.x = element_text(size=6)) +
    labs(x="Filter Stringency", y="Relative Abundance")
  
  plot3
  
  plot_grid(plot1, plot2, plot3, labels=c("A", "B", "C"))
  
  
}




