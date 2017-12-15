library(ggplot2)
library(reshape2)
library(cowplot)
setwd("~/projects/DenoiseCompare/Data_Analysis/")


Expected <- read.table("BlastResults.tsv", header = F, sep = "\t")
colnames(Expected) <- c("Study", "Pipeline", "Filter", "100% Expected", "97% Expected")

Silva <- read.table("Silva_Blast_Results.tsv", header =F, sep = "\t")
colnames(Silva) <- c("Study", "Pipeline", "Filter", "100% Database", "97% Database")


Final <- merge(Expected, Silva, all = T)

Unmatched <- read.table("Missed.txt", header = F, sep = " ")

Unmatched$nomatch <- (Unmatched[,2]-(Silva[,4] + Silva[,5]))

Final$Unmatched <- Unmatched$nomatch

Final$Percent <- ((Final$Unmatched) / (Final[,4]+Final[,5]+Final[,6]+Final[,7]+Final[,8]))*100 

LowFilt <- subset(Final, Filter=="Low")
LowFilt <- LowFilt[, -3]
LowFilt$Study <- paste(LowFilt$Study, LowFilt$Pipeline, sep = "\n")
LowFilt <- LowFilt[, -2]
LowFilt <- LowFilt[, -7]
Low_counts <- LowFilt$`100% Expected` + LowFilt$`97% Expected` + LowFilt$`100% Database` + LowFilt$`97% Database` + LowFilt$Unmatched


MedFilt <- subset(Final, Filter=="Med")

MedFilt <- MedFilt[, -3]

MedFilt$Study <- paste(MedFilt$Study, MedFilt$Pipeline, sep = "\n")
MedFilt <- MedFilt[, -2]
MedFilt <- MedFilt[, -7]
Med_counts <- MedFilt$`100% Expected` + MedFilt$`97% Expected` + MedFilt$`100% Database` + MedFilt$`97% Database` + MedFilt$Unmatched

HighFilt <- subset(Final, Filter=="High")


HighFilt <- HighFilt[, -3]
HighFilt$Study <- paste(HighFilt$Study, HighFilt$Pipeline, sep = "\n")
HighFilt <- HighFilt[, -2]
HighFilt <- HighFilt[, -7]
High_counts <- HighFilt$`100% Expected` + HighFilt$`97% Expected` + HighFilt$`100% Database` + HighFilt$`97% Database` + HighFilt$Unmatched
rbind(HighFilt$Study, MedFilt$Study, LowFilt$Study)

LowFilt_melt <- melt(data = LowFilt, id.vars = c("Study"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                     value.name = "Counts", variable.name="WorkFlow")

plotLow <- ggplot(LowFilt_melt, aes(x = Study, y = Counts, fill = WorkFlow)) +
   geom_bar(stat = "identity", show.legend = F) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
  xlab("") + scale_fill_brewer(palette="Set1")
plotLow



MedFilt_melt <- melt(data = MedFilt, id.vars = c("Study"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                     value.name = "Counts", variable.name="WorkFlow")

plotMed <- ggplot(MedFilt_melt, aes(x = Study, y = Counts, fill = WorkFlow)) +
  geom_bar(stat = "identity", show.legend = F) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
  xlab("") + scale_fill_brewer(palette="Set1")
plotMed





HighFilt_melt <- melt(data = HighFilt, id.vars = c("Study"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                     value.name = "Counts", variable.name="WorkFlow")

plotHigh <- ggplot(HighFilt_melt, aes(x = Study, y = Counts, fill = WorkFlow, labels = Study)) +
  geom_bar(stat = "identity", show.legend = F) + theme(axis.text.x = element_text(size = 5.5)) + ylab("ASV Counts") + 
  xlab("") + scale_fill_brewer(palette ="Set1")
plotHigh

plotDumby <- ggplot(HighFilt_melt, aes(x = Study, y = Counts, fill = WorkFlow)) +
  geom_bar(stat= "identity") + theme(axis.text.x = element_text(size = 5.5)) + ylab("Percent of Counts") + 
  xlab("") + scale_fill_brewer(palette = "Set1", name="Match Type")
plotDumby

Legend <- get_legend(plotDumby)
Legend
grid_plot <- plot_grid(plotHigh, plotMed, plotLow, Legend, labels=c("A", "B", "C"))
grid_plot
