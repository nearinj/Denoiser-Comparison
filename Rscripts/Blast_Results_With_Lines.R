library(ggplot2)
library(reshape2)
library(cowplot)

#Used to create stacked barcharts for the blast results

setwd("~/projects/DenoiseCompare/Data_Analysis/")
Expected <- read.table("BlastResults.tsv", header = F, sep = "\t")
colnames(Expected) <- c("Study", "Pipeline", "Filter", "100% Expected", "97% Expected")
Silva <- read.table("Silva_Blast_Results.tsv", header =F, sep = "\t")
colnames(Silva) <- c("Study", "Pipeline", "Filter", "100% Database", "97% Database")
Final <- merge(Expected, Silva, all = T)
Unmatched <- read.table("Missed.txt", header = F, sep = " ")
Unmatched$nomatch <- (Unmatched[,2]-(Silva[,4] + Silva[,5]))
Final$Unmatched <- Unmatched$nomatch

Low_Filt <- Final[Final$Filter == 'Low',]
Low_Filt <- Low_Filt[, -3]
Med_Filt <- Final[Final$Filter == 'Med',]
Med_Filt <- Med_Filt[,-3]
High_Filt <- Final[Final$Filter == 'High',]
High_Filt <- High_Filt[, -3]

PlotSamples(High_Filt, "High")
PlotSamples(Med_Filt, "Med")
PlotSamples(Low_Filt, "Low")

PlotSamples <- function(tab1, filt){
  HMP <- tab1[tab1$Study == "HMP",]  
  HMP <- HMP[,-1]
  
  Mock12 <- tab1[tab1$Study == "Mock-12",]
  Mock12 <- Mock12[,-1]
  
  Mock9 <- tab1[tab1$Study == "Mock-9",]
  Mock9 <- Mock9[,-1]
  
  Zymock <- tab1[tab1$Study == "Zymock",]
  Zymock <- Zymock[,-1]
  
  orders <- c("Unmatched", "97% Database","100% Database","97% Expected","100% Expected")
  
  Unique_expected <- read.table("Unique_counts.txt", sep="\t", header=F)
 
  
  
  HMP_melt <- melt(data = HMP, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                       value.name = "Counts", variable.name="WorkFlow")
  HMP_melt$WorkFlow <- factor(HMP_melt$WorkFlow, levels = orders)
  
  Mock12_melt <- melt(data = Mock12, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                   value.name = "Counts", variable.name="WorkFlow")
  Mock12_melt$WorkFlow <- factor(Mock12_melt$WorkFlow, levels = orders)
 
  Mock9_melt <- melt(data = Mock9, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                   value.name = "Counts", variable.name="WorkFlow")
  Mock9_melt$WorkFlow <- factor(Mock9_melt$WorkFlow, levels = orders)
  
  Zymock_melt <- melt(data = Zymock, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                   value.name = "Counts", variable.name="WorkFlow")
  Zymock_melt$WorkFlow <- factor(Zymock_melt$WorkFlow, levels = orders)
  
  
  HMP_plot <- ggplot(HMP_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[1,1], linetype="dashed")
  
  Mock12_plot <- ggplot(Mock12_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[3,1], linetype="dashed")
  
  Mock9_plot <- ggplot(Mock9_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[4,1], linetype="dashed")
  
  Zymock_plot <- ggplot(Zymock_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[2,1], linetype="dashed")
  
  grid_plot <- plot_grid(HMP_plot, Mock12_plot, Mock9_plot, Zymock_plot, labels=c("A", "B", "C", "D"))
  save_plot(paste("Blast_Comp_Figs/ASV Comparison_", filt, ".png", sep=""), grid_plot, base_aspect_ratio = 2)
}
