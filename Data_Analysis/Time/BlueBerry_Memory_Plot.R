#used to plot memory usage foor each pipeline

library(ggplot2)

setwd("~/GitHub_Repos/Denoiser Comparison/Data_Analysis/Time/")

MasterTable <- read.table("Blueberry_time_data.tsv", sep="\t", header=T)

#get Deblur memory
deblur_memory <- MasterTable[grep("Deblur", MasterTable$Pipe),c("Filt", "TimeCMD","Memory.kb.")]
deblur_memory$Pipe <- "Deblur"
#Remove 40k reads as the read depth was deep enough to support this rarefaction
deblur_memory <- deblur_memory[-5,]


#get DADA2 memory

MemoryDataInference <- MasterTable[which(MasterTable$TimeCMD == "timeI"),]
MemoryDataChimera <- MasterTable[which(MasterTable$TimeCMD == "timeCT"),]
rownames(MemoryDataInference) <- paste(MemoryDataInference$Filt, MemoryDataInference$Pipe, sep="_")
rownames(MemoryDataChimera) <- paste(MemoryDataChimera$Filt, MemoryDataChimera$Pipe, sep="_")
dada2_memory <- MemoryDataInference

dada2_memory <- dada2_memory[,-c(4,5)]
dada2_memory <- dada2_memory[-5,]
dada2_memory$Pipe <- "DADA2"

#get UNOISE memory


MemoryDataD <- MasterTable[which(MasterTable$TimeCMD == "Dtime"),]
MemoryDataUni <- MasterTable[which(MasterTable$TimeCMD == "Unitime"),]
rownames(MemoryDataD) <- paste(MemoryDataD$Filt, MemoryDataD$Pipe, sep="_")
rownames(MemoryDataUni) <- paste(MemoryDataUni$Filt, MemoryDataUni$Pipe, sep="_")
Unoise_Memory <- MemoryDataUni
Unoise_Memory <- Unoise_Memory[,-c(4,5)]
Unoise_Memory <- Unoise_Memory[-5,]
Unoise_Memory$Pipe <- "UNOISE3"

#Make one large table
recombined_Memory <- rbind(deblur_memory, dada2_memory, Unoise_Memory)
colnames(recombined_Memory)[4] <- "Pipeline"

#set the bar order
orders <- c("UNOISE3 325k", "UNOISE3 648k", "UNOISE3 1287k", "UNOISE3 1926k", 
            "Deblur 325k", "Deblur 648k", "Deblur 1287k", "Deblur 1926k",
            "DADA2 325k", "DADA2 648k", "DADA2 1287k", "DADA2 1926k")
#rename column to be used as labels in graph
recombined_Memory$Filt <- paste(recombined_Memory$Pipe, c("325k", "648k", "1287k", "1926k"), sep = " ")
group.colors <- c(Deblur = "#333BFF", DADA2 = "#CC6600", UNOISE3 = "#9633FF")
recombined_Memory$Filt <- factor(recombined_Memory$Filt, levels = orders)

#make the plot
Memoryplot <- ggplot(recombined_Memory, aes(x= Filt, y=Memory.kb., fill=Pipeline)) + geom_bar(stat = "identity") + scale_fill_manual(values=group.colors) +
  xlab("") + theme(axis.text = element_text(size=10)) + ylab("Memory used in Kb") + scale_y_continuous(breaks=c(1000000,2000000,3000000,4000000,5000000))
Memoryplot
