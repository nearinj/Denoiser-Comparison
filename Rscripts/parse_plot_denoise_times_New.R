library(ggplot2)


#Read in time table
timeData <- read.table("~/projects/DenoiseCompare/Data_Analysis/time/Time_Table.tsv", 
                       header = TRUE, 
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       )
#Get Deblur Times
deblur_times <- timeData[grep("Deblur", timeData$Pipe),c("Filt", "TimeCMD", "Real", "User", "Sys")]
deblur_times$Pipe <- "Deblur"
deblur_times <- deblur_times[-5,]


#Combine Dada2 Times

timeDataInference <- timeData[which(timeData$TimeCMD == "timeI"),]
timeDataChimera <- timeData[which(timeData$TimeCMD == "timeCT"),]

rownames(timeDataInference) <- paste(timeDataInference$Filt, timeDataInference$Pipe, sep="_")
rownames(timeDataChimera) <- paste(timeDataChimera$Filt, timeDataChimera$Pipe, sep="_")

dada2_time <- timeDataInference
dada2_time$User <- timeDataInference$User + timeDataChimera[rownames(timeDataInference), "User"]
dada2_time <- dada2_time[-5,]

# Combine Unoise Times

timeDataD <- timeData[which(timeData$TimeCMD == "Dtime"),]
timeDataUni <- timeData[which(timeData$TimeCMD == "Unitime"),]

rownames(timeDataD) <- paste(timeDataD$Filt, timeDataD$Pipe, sep="_")
rownames(timeDataUni) <- paste(timeDataUni$Filt, timeDataUni$Pipe, sep="_")

Unoise_time <- timeDataD
Unoise_time$User <- timeDataD$User + timeDataUni[rownames(timeDataD), "User"]
Unoise_time <- Unoise_time[-5, ]

recombined_time <- rbind(deblur_times, dada2_time, Unoise_time)
orders <- c("Unoise 515k", "Unoise 1023k", "Unoise 1939k", "Unoise 2495k", 
                                  "Deblur 515k", "Deblur 1023k", "Deblur 1939k", "Deblur 2495k",
                                  "Dada 515k", "Dada 1023k", "Dada 1939k", "Dada 2495k")
colnames(recombined_time)[6]<- "Pipeline"
recombined_time$Filt <- paste(recombined_time$Pipeline, c("515k", "1023k", "1939k", "2495k"), sep = " ")


#barplot(log10(recombined_time$User), names.arg = argnames, las=2, ylab = "log10(User time in s)", main = "Speed of the Different Pipelines", )

group.colors <- c(Deblur = "#333BFF", Dada = "#CC6600", Unoise = "#9633FF")

recombined_time$Filt <- factor(recombined_time$Filt, levels = orders)

timeplot <- ggplot(recombined_time, aes(x= Filt, y=log10(User), fill=Pipeline)) + geom_bar(stat = "identity") + scale_fill_manual(values=group.colors) +
  xlab("") + theme(axis.text = element_text(size=10)) + ylab("User Time in Seconds") + scale_y_continuous(breaks=c(0,1,2,3,4,5), labels=c("0","10","100","1000","10000","100000"))
timeplot
