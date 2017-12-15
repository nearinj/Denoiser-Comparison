library(dada2)

setwd("~/projects/DenoiseCompare_Out/Blueberry/med/deblurP/final/")
seqtab <- readRDS("~/projects/DenoiseCompare_Out/Blueberry/med/dada2/seqtab_final.rds")


deblurtable <- read.table("table.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
deblurID <- read.table("seqID.tsv", header=FALSE, sep="\t")

target <- deblurtable$OTU.ID

deblurIDtemp <- deblurID[match(target, deblurID$V1),]
dim(deblurID)

deblurtable$OTU.ID <- deblurIDtemp$V2
deblurtablelong <- t(deblurtable)


colnames(deblurtablelong) = deblurtablelong[1,]
deblurtablelong <- deblurtablelong[-1,]
dim(deblurtablelong)
class(deblurtablelong) <- "numeric"

taxa <- assignTaxonomy(deblurtablelong, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = T)

saveRDS(taxa, file="tax_final.rds")
saveRDS(deblurtablelong, file="seqtab_final.rds")

setwd("../../Unoise/")
unoisetable <- read.table("otutab.txt", header=TRUE, sep="\t")
unoiseID <- read.table("otuID.tsv", header=F, sep="\t")

target <- unoisetable$OTU.ID
unoiseIDtemp <- unoiseID[match(target, unoiseID$V1),]
unoisetable$OTU.ID <- unoiseIDtemp$V2

rownames(unoisetable) = unoisetable[,1]
unoisetable <- unoisetable[,-1]
dim(unoisetable)
colnames(unoisetable) <- gsub("_*","",colnames(unoisetable))
colnames(unoisetable) <- gsub("R1filt","",colnames(unoisetable))

unoiselong <- t(unoisetable)
class(unoiselong) <- "numeric"
str(unoiselong)


taxa <- assignTaxonomy(unoiselong, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = T)
saveRDS(taxa, file="tax_final.rds")
saveRDS(unoiselong, file="seqtab_final.rds")
