#remove_low blast results chart for HMP... do Zymock and Fungal as well..
library(ggplot2)
library(cowplot)
library(reshape2)


setwd("/home/jacob/projects/DenoiseCompare/Data_Analysis/")

# Read in results of blasting study seqs to expected sequences after removal of low confidence seqs.
Expected <- read.table("remove_low_BlastResults.tsv", header = F, sep = "\t")
colnames(Expected) <- c("Study", "Pipeline", "Filter", "100% Expected", "97% Expected")

# Read in results of blasting study seqs to SILVA database after removal of low confidence seqs.
Silva <- read.table("remove_low_Silva_Blast_Results.tsv", header =F, sep = "\t")
colnames(Silva) <- c("Study", "Pipeline", "Filter", "100% Database", "97% Database")

#merge the two tables together
Final <- merge(Expected, Silva, all = T)

FinalHMP <- subset(Final, Final$Study == "HMP")
FinalZymock <- subset(Final, Final$Study == "Zymock")
FinalMock9 <- subset(Final, Final$Study == "Mock9")

FinalHMP$Pipeline <- gsub("Dada", "DADA2", FinalHMP$Pipeline)
FinalHMP$Pipeline <- gsub("Unoise", "UNOISE3", FinalHMP$Pipeline)
#make the unmatched column
FinalHMP$Unmatched <- c(0, 0, 4)

FinalHMP_melt <- melt(FinalHMP, measure.vars=c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"), id.vars="Pipeline", value.name="Counts", variable.name="Workflow")

Unique_expected <-read.table("Unique_counts.txt", sep="\t", header=F)

orders <- c("Unmatched", "97% Database","100% Database","97% Expected","100% Expected")

FinalHMP_melt$Workflow <- factor(FinalHMP_melt$Workflow, levels = orders)
HMP_plot <- ggplot(FinalHMP_melt, aes(x=Pipeline, y=Counts, fill= Workflow))+ 
  geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 14)) + ylab("ASV Counts") + 
  xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[1,1], linetype="dashed")
HMP_plot


FinalZymock$Pipeline <- gsub("Dada", "DADA2", FinalZymock$Pipeline)
FinalZymock$Pipeline <- gsub("Unoise", "UNOISE3", FinalZymock$Pipeline)
#make the unmatched column Dada, Deblur, Unoise, 19-1
FinalZymock$Unmatched <- c(0, 0, 18)

FinalZymock_melt <- melt(FinalZymock, measure.vars=c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"), id.vars="Pipeline", value.name="Counts", variable.name="Workflow")


FinalZymock_melt$Workflow <- factor(FinalZymock_melt$Workflow, levels = orders)
Zymock_plot <- ggplot(FinalZymock_melt, aes(x=Pipeline, y=Counts, fill= Workflow))+ 
  geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 14)) + ylab("ASV Counts") + 
  xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[2,1], linetype="dashed")
Zymock_plot

FinalMock9$Pipeline <- gsub("Dada", "DADA2", FinalMock9$Pipeline)
FinalMock9$Pipeline <- gsub("Unoise", "UNOISE3", FinalMock9$Pipeline)
#make the unmatched column Dada 7-5, Deblur 16-10, Unoise, 25-15
FinalMock9$Unmatched <- c(2, 6, 10)

FinalMock9_melt <- melt(FinalMock9, measure.vars=c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"), id.vars="Pipeline", value.name="Counts", variable.name="Workflow")


FinalMock9_melt$Workflow <- factor(FinalMock9_melt$Workflow, levels = orders)
Mock9_plot <- ggplot(FinalMock9_melt, aes(x=Pipeline, y=Counts, fill= Workflow))+ 
  geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 14)) + ylab("ASV Counts") + 
  xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[2,1], linetype="dashed")
Mock9_plot


grid_plot <- plot_grid(HMP_plot, Mock9_plot, Zymock_plot, labels=c("A", "B", "C"))
grid_plot



