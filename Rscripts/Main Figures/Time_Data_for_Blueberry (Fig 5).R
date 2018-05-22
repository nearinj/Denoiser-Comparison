library(ggplot2)

setwd("~/GitHub_Repos/Denoiser Comparison/Data_Analysis/Time/")

MasterTable <- read.table("Blueberry_time_data.tsv", sep="\t", header=T)

#get Deblur times
deblur_times <- MasterTable[grep("Deblur", MasterTable$Pipe),c("Filt", "TimeCMD","User", "Sys")]
deblur_times$Pipe <- "Deblur"
deblur_times <- deblur_times[,-4]
deblur_times <- deblur_times[-5,]

#Combined Dada2
timeDataInference <- MasterTable[which(MasterTable$TimeCMD == "timeI"),]
timeDataChimera <- MasterTable[which(MasterTable$TimeCMD == "timeCT"),]
rownames(timeDataInference) <- paste(timeDataInference$Filt, timeDataInference$Pipe, sep="_")
rownames(timeDataChimera) <- paste(timeDataChimera$Filt, timeDataChimera$Pipe, sep="_")
dada2_time <- timeDataInference
dada2_time$User <- timeDataInference$User + timeDataChimera[rownames(timeDataInference), "User"]
dada2_time <- dada2_time[,-c(5,6)]
dada2_time <- dada2_time[-5,]
dada2_time$Pipe <- "DADA2"
#Combined UnoiseTimes

timeDataD <- MasterTable[which(MasterTable$TimeCMD == "Dtime"),]
timeDataUni <- MasterTable[which(MasterTable$TimeCMD == "Unitime"),]
rownames(timeDataD) <- paste(timeDataD$Filt, timeDataD$Pipe, sep="_")
rownames(timeDataUni) <- paste(timeDataUni$Filt, timeDataUni$Pipe, sep="_")
Unoise_time <- timeDataD
Unoise_time$User <- timeDataD$User + timeDataUni[rownames(timeDataD), "User"]
Unoise_time <- Unoise_time[,-c(5,6)]
Unoise_time <- Unoise_time[-5,]
Unoise_time$Pipe <- "UNOISE3"

recombined_time <- rbind(deblur_times, dada2_time, Unoise_time)
colnames(recombined_time)[4] <- "Pipeline"

orders <- c("UNOISE3 325k", "UNOISE3 648k", "UNOISE3 1287k", "UNOISE3 1926k", 
            "Deblur 325k", "Deblur 648k", "Deblur 1287k", "Deblur 1926k",
            "DADA2 325k", "DADA2 648k", "DADA2 1287k", "DADA2 1926k")



recombined_time$Filt <- paste(recombined_time$Pipe, c("325k", "648k", "1287k", "1926k"), sep = " ")


#barplot(log10(recombined_time$User), names.arg = argnames, las=2, ylab = "log10(User time in s)", main = "Speed of the Different Pipelines", )

group.colors <- c(Deblur = "#333BFF", DADA2 = "#CC6600", UNOISE3 = "#9633FF")

recombined_time$Filt <- factor(recombined_time$Filt, levels = orders)

timeplot <- ggplot(recombined_time, aes(x= Filt, y=log10(User), fill=Pipeline)) + geom_bar(stat = "identity") + scale_fill_manual(values=group.colors) +
  xlab("") + theme(axis.text = element_text(size=10)) + ylab("User Time in Seconds") + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("0","10","100","1000","10000","100000"))
timeplot


#grid_plot <- plot_grid(timeplot, Memoryplot, ncol = 1, labels = c("A", "B"))
#grid_plot
