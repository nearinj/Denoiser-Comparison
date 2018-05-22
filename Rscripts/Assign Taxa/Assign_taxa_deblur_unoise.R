#used to assign taxa to the blueberry data generated from the deblur pipeline
library(dada2)

setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/deblurP/final/")

deblurtable <- read.table("table.tsv", 
                          header=TRUE, 
                          sep="\t", 
                          stringsAsFactors = F, 
                          skip=1, 
                          comment.char="",
                          row.names=1)

#generated from rename_seqs.py

deblurID <- read.table("seqID.tsv", header=FALSE, sep="\t")

target <- rownames(deblurtable)

deblurIDtemp <- deblurID[match(target, deblurID$V1),]
dim(deblurID)

rownames(deblurtable) <- deblurIDtemp$V2
deblurtablelong <- t(deblurtable)

dim(deblurtablelong)

deblur_taxa <- assignTaxonomy(deblurtablelong, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 40)

saveRDS(deblur_taxa, file="tax_final.rds")
saveRDS(deblurtablelong, file="seqtab_final.rds")

### Assign taxonomy to UNOISE sequences.
setwd("/home/jacob/projects/DenoiseCompare_Out/Blueberry/med/Unoise/")
unoisetable <- read.table("otutab.txt",                           
                          header=TRUE, 
                          sep="\t", 
                          stringsAsFactors = F, 
                          comment.char="",
                          row.names=1)

unoiseID <- read.table("otuID.tsv", header=F, sep="\t")

target <- rownames(unoisetable)

unoiseIDtemp <- unoiseID[match(target, unoiseID$V1),]
rownames(unoisetable) <- unoiseIDtemp$V2
colnames(unoisetable) <- gsub("_R1_filt$","",colnames(unoisetable))
unoiselong <- t(unoisetable)

unoise_taxa <- assignTaxonomy(unoiselong, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 40)
saveRDS(unoise_taxa, file="tax_final.rds")
saveRDS(unoiselong, file="seqtab_final.rds")
