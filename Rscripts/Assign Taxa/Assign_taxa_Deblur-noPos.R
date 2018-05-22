#assign taxa to deblur without postive filtering


library(dada2)


##################################### biscuit data

setwd("projects/DenoiseCompare_Out/biscuit/med/Deblur_NoPos/output/")

biscuit_table <- read.table("table.tsv", header=T, sep="\t", comment.char="", skip=1, row.names=1)
colnames(biscuit_table) <- gsub("_.", "", colnames(biscuit_table))

biscuit_table_long <- t(biscuit_table)
dim(biscuit_table_long)

biscuit_taxa <- assignTaxonomy(biscuit_table_long, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 20)
dim(biscuit_taxa)


saveRDS(biscuit_taxa, file="tax_final.rds")
saveRDS(biscuit_table_long, file="seqtab_final.rds")


############################################# Exercise data

setwd("~/projects/DenoiseCompare_Out/Exercise_med/med/Deblur_NoPos/output/")

exercise_table <- read.table("table.tsv", header=T, sep="\t", comment.char="", skip=1, row.names=1)
colnames(exercise_table) <- gsub("_.", "", colnames(exercise_table))

exercise_table_long <- t(exercise_table)
dim(exercise_table_long)

exercise_taxa <- assignTaxonomy(exercise_table_long, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 20)
dim(exercise_taxa)


saveRDS(exercise_taxa, file="tax_final.rds")
saveRDS(exercise_table_long, file="seqtab_final.rds")


############################################ Blueberry data

setwd("~/projects/DenoiseCompare_Out/Blueberry/med/Deblur_NoPos/output/")

blueberry_table <- read.table("table.tsv", header=T, sep="\t", comment.char="", skip=1, row.names=1)
colnames(blueberry_table) <- gsub("_.", "", colnames(blueberry_table))

blueberry_table_long <- t(blueberry_table)
dim(blueberry_table_long)

blueberry_taxa <- assignTaxonomy(blueberry_table_long, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 20)
dim(blueberry_taxa)


saveRDS(blueberry_taxa, file="tax_final.rds")
saveRDS(blueberry_table_long, file="seqtab_final.rds")


### once this is done convert bioms to tsv, add headers to the sequence names convert back.


