#combined tsv into single tsv to be converted into a biom file for further analysis

#for Exercise

setwd("~/projects/DenoiseCompare_Out/Exercise_med/med/COMBINED/biom/")

Dada2_table <- read.table("Dada2.tsv", 
                          header=TRUE, 
                          sep="\t", 
                          stringsAsFactors = F, 
                          skip=1, 
                          comment.char=""
                          )
Deblur_table <- read.table("Deblur.tsv", 
                          header=TRUE, 
                          sep="\t", 
                          stringsAsFactors = F, 
                          skip=1, 
                          comment.char=""
                          )
Unoise_table <- read.table("Unoise.tsv", 
                          header=TRUE, 
                          sep="\t", 
                          stringsAsFactors = F, 
                          skip=1, 
                          comment.char=""
                          )

#make sure they all contain the same number of samples
Dada2_table <- Dada2_table[, colnames(Deblur_table)]
Unoise_table <- Unoise_table[, colnames(Deblur_table)]


#add prefix to seqs
Dada2_table$X.OTU.ID <- paste("dada", Dada2_table$X.OTU.ID, sep="_")
Deblur_table$X.OTU.ID <- paste("deblur", Deblur_table$X.OTU.ID, sep="_")
Unoise_table$X.OTU.ID <- paste("unoise", Unoise_table$X.OTU.ID, sep="_")


#add prefix to sample
colnames(Dada2_table)[-1] <- paste("Dada", colnames(Dada2_table)[-1], sep="_")
colnames(Deblur_table)[-1] <- paste("Deblur", colnames(Deblur_table)[-1], sep="_")
colnames(Unoise_table)[-1] <- paste("Unoise_", colnames(Unoise_table)[-1], sep="_")


#combined tables together
Combined <- merge(Dada2_table, Deblur_table, all=T)
Combined <- merge(Combined, Unoise_table, all=T)

#fill in Na as 0
Combined[is.na(Combined)] <- 0
#rename otu table
colnames(Combined)[1] <- "OTU.ID"

#write TSV file
write.table(Combined, file="Combined_Exercise.tsv", sep="\t", quote=F, col.names = NA)


#write the taxonomy files to be added to the biom table generated from the above table


#read taxa files in

Deblur_taxa <- read.table("../../deblurP/final/taxa_metadata.txt", sep="\t", header=T, comment.char="")
Unoise_taxa <- read.table("../../Unoise/taxa_metadata.txt", sep="\t", header=T, comment.char="")
Dada2_taxa <- read.table("../../dada2/taxa_metadata.txt", sep="\t", header=T, comment.char="")


#add prefixs to seq names
Deblur_taxa$X.OTU.ID <- paste("deblur", Deblur_taxa$X.OTU.ID, sep="_")
Unoise_taxa$X.OTU.ID <- paste("unoise", Unoise_taxa$X.OTU.ID, sep="_")
Dada2_taxa$X.OTU.ID <- paste("dada", Dada2_taxa$X.OTU.ID, sep="_")


#merge into one large taxa file
taxa <- merge(Deblur_taxa, Unoise_taxa, all=T)
taxa <- merge(taxa, Dada2_taxa, all=T)

#write to a tsv file
write.table(taxa, file="taxa_Combined_Exercise.tsv", sep="\t", quote=F, col.names = NA)
