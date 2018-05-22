library(ggplot2)
library(reshape2)
library(cowplot)

#Used to create stacked barcharts for the blast results

setwd("~/projects/DenoiseCompare/Data_Analysis/fixed_ASV_analysis/")

# Read in results of blasting study seqs to expected sequences.
Expected <- read.table("BlastResults.tsv", header = F, sep = "\t", stringsAsFactors = F)
colnames(Expected) <- c("Study", "Pipeline", "Filter", "100% Expected", "97% Expected")

# Read in results of blasting study seqs to SILVA database.
Silva <- read.table("Silva_Blast_Results.tsv", header =F, sep = "\t", stringsAsFactors = F)
colnames(Silva) <- c("Study", "Pipeline", "Filter", "100% Database", "97% Database", "Unmatched")

Final <- merge(Expected, Silva, all = T)
#set rownames for later use
rownames(Final) <- paste(Final$Study, Final$Pipeline, Final$Filter, sep="_")

#get unmatched amount

total_query <- read.table("Total_query.tsv", sep="\t")
colnames(total_query) <- c("Study", "Pipeline", "Filter", "Total_Querys")
rownames(total_query) <- paste(total_query$Study, total_query$Pipeline, total_query$Filter, sep="_")
#set total_query to be same order as Final
total_query <- total_query[rownames(Final),]
total_hit <- rowSums(Final[,c(4,5,6,7,8)])
total_hit == total_query$Total_Querys
#all sequences were hit!



# Subset final table to different filter settings.
Low_Filt <- Final[Final$Filter == 'Low',]
Low_Filt <- Low_Filt[, -3]
Med_Filt <- Final[Final$Filter == 'Med',]
Med_Filt <- Med_Filt[,-3]
High_Filt <- Final[Final$Filter == 'High',]
High_Filt <- High_Filt[, -3]

#Subset med to orignial and with open-ref

Med_Filt_orig <- Med_Filt[-(which(Med_Filt$Pipeline=="Open")),]


# Plot samples (need to load in below function first)
PlotSamples(High_Filt, "High")
plot_with_open(Med_Filt, "Med")
PlotSamples(Low_Filt, "Low")

#plot Orig
PlotSamples(Med_Filt_orig, "Med")

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
  
  #grid_plot <- 
  plot_grid(HMP_plot, Mock12_plot, Mock9_plot, Zymock_plot, labels=c("A", "B", "C", "D"))
  
  #save_plot(paste("Blast_Comp_Figs/ASV Comparison_", filt, ".png", sep=""), grid_plot, base_aspect_ratio = 2)
  
}


# plot with Open

library(grid)

plot_with_open <- function(tab1, filt){
  
  HMP <- tab1[tab1$Study == "HMP",]  
  HMP <- HMP[,-1]
  HMP_Open <- HMP[which(HMP$Pipeline=="Open"),]
  HMP <- HMP[-which(HMP$Pipeline=="Open"),]
  
  Mock12 <- tab1[tab1$Study == "Mock-12",]
  Mock12 <- Mock12[,-1]
  Mock12_Open <- Mock12[which(Mock12$Pipeline=="Open"),]
  Mock12 <- Mock12[-which(Mock12$Pipeline=="Open"),]
  
  Mock9 <- tab1[tab1$Study == "Mock-9",]
  Mock9 <- Mock9[,-1]
  Mock9_Open <- Mock9[which(Mock9$Pipeline == "Open"),]
  Mock9 <- Mock9[-which(Mock9$Pipeline=="Open"),]
  
  Zymock <- tab1[tab1$Study == "Zymock",]
  Zymock <- Zymock[,-1]
  Zymock_Open <- Zymock[which(Zymock$Pipeline == "Open"),]
  Zymock <- Zymock[-which(Zymock$Pipeline == "Open"),]
  
  orders <- c("Unmatched", "97% Database","100% Database","97% Expected","100% Expected")
  
  Unique_expected <- read.table("Unique_counts.txt", sep="\t", header=F)
  
  HMP_melt <- melt(data = HMP, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                   value.name = "Counts", variable.name="WorkFlow")
  HMP_melt$WorkFlow <- factor(HMP_melt$WorkFlow, levels = orders)
  ########
  HMP_Open_melt <- melt(data = HMP_Open, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                        value.name = "Counts", variable.name="WorkFlow")
  HMP_Open_melt$WorkFlow <- factor(HMP_Open_melt$WorkFlow, levels=orders) 
  #########
  Mock12_melt <- melt(data = Mock12, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                      value.name = "Counts", variable.name="WorkFlow")
  Mock12_melt$WorkFlow <- factor(Mock12_melt$WorkFlow, levels = orders)
  #########
  Mock12_Open_melt <- melt(data = Mock12_Open, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                           value.name = "Counts", variable.name="WorkFlow")
  Mock12_Open_melt$WorkFlow <- factor(Mock12_Open_melt$WorkFlow, levels = orders)
  ##########
  Mock9_melt <- melt(data = Mock9, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                     value.name = "Counts", variable.name="WorkFlow")
  Mock9_melt$WorkFlow <- factor(Mock9_melt$WorkFlow, levels = orders)
  ###########
  Mock9_Open_melt <- melt(data = Mock9_Open, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                          value.name = "Counts", variable.name="WorkFlow")
  Mock9_Open_melt$WorkFlow <- factor(Mock9_Open_melt$WorkFlow, levels = orders)
  ###################
  Zymock_melt <- melt(data = Zymock, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                      value.name = "Counts", variable.name="WorkFlow")
  Zymock_melt$WorkFlow <- factor(Zymock_melt$WorkFlow, levels = orders)
  ################
  Zymock_Open_melt <- melt(data = Zymock_Open, id.vars = c("Pipeline"), measured.vars = c("100% Expected", "97% Expected", "100% Database", "97% Database", "Unmatched"),
                           value.name = "Counts", variable.name="WorkFlow")
  Zymock_Open_melt$WorkFlow <- factor(Zymock_Open_melt$WorkFlow, levels = orders)
  ##############################
  
  ################ HMP PLOT
  HMP_plot <- ggplot(HMP_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12), legend.position = "none") + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[1,1], linetype="dashed") + ylim(0,65)
  ### Open plot
  HMP_open_plot <- ggplot(HMP_Open_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12)) + ylab("") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[1,1], linetype="dashed") + ylim(0,65)
  ### Combine them together
  Final_HMP_plot <- plot_grid(HMP_plot, HMP_open_plot, rel_widths = (c(1.5,1)), align="hv")
  
  ############### Mock12 PLOT
  Mock12_plot <- ggplot(Mock12_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12), legend.position = "none") + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[3,1], linetype="dashed") 
  
  Mock12_Open_plot <- ggplot(Mock12_Open_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12)) + ylab("") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[3,1], linetype="dashed")
  
  Final_Mock12_plot <- plot_grid(Mock12_plot, Mock12_Open_plot, rel_widths = c(1.5,1), align="hv")
  
  
  #################### Mock9 PLOT
  Mock9_plot <- ggplot(Mock9_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12), legend.position = "none") + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[4,1], linetype="dashed") + ylim(0,45)
  
  Mock9_open_plot <- ggplot(Mock9_Open_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12)) + ylab("") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[4,1], linetype="dashed") + ylim(0,45)
  
  Final_Mock9_plot <- plot_grid(Mock9_plot, Mock9_open_plot, rel_widths = c(1.5,1), align="hv")
  
  Zymock_plot <- ggplot(Zymock_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12), legend.position = "none") + ylab("ASV Counts") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[2,1], linetype="dashed") + ylim(0,45)
  
  Zymock_Open_plot <- ggplot(Zymock_Open_melt, aes(x = Pipeline, y = Counts, fill = WorkFlow)) +
    geom_bar(stat = "identity", show.legend = T) + theme(axis.text.x = element_text(size = 12)) + ylab("") + 
    xlab("") + scale_fill_brewer(palette="Set1") + geom_hline(yintercept = Unique_expected[2,1], linetype="dashed") + ylim(0,45)
  
  Final_Zymock_plot <- plot_grid(Zymock_plot, Zymock_Open_plot, rel_widths = c(1.5,1), align="hv")
  
  #grid_plot <- 
  plot_grid(Final_HMP_plot, Final_Mock12_plot, Final_Mock9_plot, Final_Zymock_plot, labels=c("A", "B", "C", "D"))
  
}







### Filtered ASV figures
setwd("~/projects/DenoiseCompare/Data_Analysis/fixed_ASV_analysis/remove_low/")

# Read in results of blasting study seqs to expected sequences.
R_Expected <- read.table("BlastResults.tsv", header = F, sep = "\t", stringsAsFactors = F)
colnames(R_Expected) <- c("Study", "Pipeline", "Filter", "100% Expected", "97% Expected")

# Read in results of blasting study seqs to SILVA database.
R_Silva <- read.table("Silva_Blast_Results.tsv", header =F, sep = "\t", stringsAsFactors = F)
colnames(R_Silva) <- c("Study", "Pipeline", "Filter", "100% Database", "97% Database", "Unmatched")

R_Final <- merge(R_Expected, R_Silva, all = T)
#set rownames for later use
rownames(R_Final) <- paste(R_Final$Study, R_Final$Pipeline, R_Final$Filter, sep="_")

R_Final$Study <- gsub("Mock9", "Mock-9", R_Final$Study)
R_Final$Study <- gsub("Mock12", "Mock-12", R_Final$Study)
#get unmatched amount

R_total_query <- read.table("Total_query.tsv", sep="\t")
colnames(R_total_query) <- c("Study", "Pipeline", "Filter", "Total_Querys")
rownames(R_total_query) <- paste(R_total_query$Study, R_total_query$Pipeline, R_total_query$Filter, sep="_")
#set total_query to be same order as Final
R_total_query <- R_total_query[rownames(R_Final),]
R_total_hit <- rowSums(R_Final[,c(4,5,6,7,8)])
R_total_hit == R_total_query$Total_Querys
#all sequences were hit!






# Subset final table to different filter settings.

R_Med_Filt <- R_Final[R_Final$Filter == 'Med',]
R_Med_Filt <- R_Med_Filt[,-3]


#Subset med to orignial and with open-ref

R_Med_Filt_orig <- R_Med_Filt[-(which(R_Med_Filt$Pipeline=="Open")),]

# Plot samples (need to load in below function first)

plot_with_open(R_Med_Filt, "Med")



#plot Orig
PlotSamples(Med_Filt_orig, "Med")
#save plot as 11x19.5

#