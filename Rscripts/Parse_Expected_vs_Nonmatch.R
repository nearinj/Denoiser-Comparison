BlastResults <- read.table("~/projects/DenoiseCompare/Data_Analysis/BlastResults.tsv", 
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE)
ColumnNames <- c("Sample", "Pipeline", "Filter", "Expected 100%", "Expected 97%")
colnames(BlastResults) <- ColumnNames
SilvaResults <- read.table("~/projects/DenoiseCompare/Data_Analysis/Silva_Blast_Results.tsv",
                              header = FALSE,
                              sep = "\t")
SilvaNames <- c("Sample", "Pipeline", "Filter", "Silva 100%", "Silva 97%")
colnames(SilvaResults) <- SilvaNames

FinalResults <- merge(BlastResults, SilvaResults, by=c("Sample", "Pipeline", "Filter"))


UnmatchedTotal <- SilvaResults$`Silva 100%`+ SilvaResults$`Silva 97%`

Total_nonmatch <- read.table("~/projects/DenoiseCompare/Data_Analysis/Missed.txt",
                             header = FALSE,
                             sep = " ")
Missed <- Total_nonmatch$V2-UnmatchedTotal

FinalResults$No_Match <- Missed

TotalSequences <- FinalResults$`Expected 100%`+FinalResults$`Expected 97%`+FinalResults$`Silva 100%`+FinalResults$`Silva 97%`+FinalResults$No_Match
FinalResults$Total_Sequences <- TotalSequences
