#Assign Taxa for all Open_Ref OTU from real datasets.



library(dada2)


setwd("~/projects/DenoiseCompare_Out/biscuit/med/OpenRef_OTU/remove_one/final/")


biscuit_table <- read.table("table.tsv", 
                            header=TRUE, 
                            sep="\t", 
                            stringsAsFactors = F, 
                            skip=1, 
                            comment.char="",
                            row.names=1)
biscuit_ID <- read.table("seqID.tsv", header=F, sep="\t")

target <- rownames(biscuit_table)

biscuitIDtemp <- biscuit_ID[match(target,biscuit_ID$V1),]


#sanity checks
#check to make sure row names for biscuit table are the same as the replacement.
#check if any don't equal each other
length(which(rownames(biscuit_table) == biscuitIDtemp$V1))

rownames(biscuit_table) <- biscuitIDtemp$V2



#sanity checks
rownames(biscuit_ID) <- biscuit_ID$V1
rownames(biscuitIDtemp) <- biscuitIDtemp$V1
#sanity check to make sure that nothing was switched around.
biscuit_ID["00d869123bb648c63108abf9c68c0fbf3c415a18", "V2"] == biscuitIDtemp["00d869123bb648c63108abf9c68c0fbf3c415a18", "V2"]
biscuit_ID["12574", "V2"] == biscuitIDtemp["12574", "V2"]
#looks good.

biscuit_table_long <- t(biscuit_table)
dim(biscuit_table_long)
biscuit_taxa <- assignTaxonomy(biscuit_table_long, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 20)
dim(biscuit_taxa)

saveRDS(biscuit_taxa, file="tax_final.rds")
saveRDS(biscuit_table_long, file="seqtab_final.rds")



#Blueberry Data set


setwd("~/projects/DenoiseCompare_Out/Blueberry/med/OpenRef_OTU/remove_one/final/")


blueberry_table <-  read.table("table.tsv", 
                               header=TRUE, 
                               sep="\t", 
                               stringsAsFactors = F, 
                               skip=1, 
                               comment.char="",
                               row.names=1)
blueberry_ID <- read.table("seqID.tsv", header=F, sep="\t")

target <- rownames(blueberry_table)

blueberryIDtemp <- blueberry_ID[match(target,blueberry_ID$V1),]


#sanity checks
#check to make sure row names for biscuit table are the same as the replacement.
#check if any don't equal each other
length(which(rownames(blueberry_table) == blueberryIDtemp$V1))

rownames(blueberry_table) <- blueberryIDtemp$V2

length(which(rownames(blueberry_table) == blueberryIDtemp$V2))


#sanity checks
rownames(blueberry_ID) <- blueberry_ID$V1
rownames(blueberryIDtemp) <- blueberryIDtemp$V1
#sanity check to make sure that nothing was switched around.
blueberry_ID["0000b13b5090bff64ce9c9b8d74a5a770349912d", "V2"] == blueberryIDtemp["0000b13b5090bff64ce9c9b8d74a5a770349912d", "V2"]


#looks good.

blueberry_table_long <- t(blueberry_table)
dim(blueberry_table_long)
blueberry_taxa <- assignTaxonomy(blueberry_table_long, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 20)
dim(blueberry_taxa)

saveRDS(blueberry_taxa, file="tax_final.rds")
saveRDS(blueberry_table_long, file="seqtab_final.rds")



#taxonomy for Exercise dataset

setwd("~/projects/DenoiseCompare_Out/Exercise_med/med/OpenRef_OTU/remove_one/final/")


Exercise_table <-  read.table("table.tsv", 
                               header=TRUE, 
                               sep="\t", 
                               stringsAsFactors = F, 
                               skip=1, 
                               comment.char="",
                               row.names=1)
Exercise_ID <- read.table("seqID.tsv", header=F, sep="\t")

target <- rownames(Exercise_table)

ExerciseIDtemp <- Exercise_ID[match(target,Exercise_ID$V1),]


#sanity checks
#check to make sure row names for biscuit table are the same as the replacement.
#check if any don't equal each other
length(which(rownames(Exercise_table) == ExerciseIDtemp$V1))

rownames(Exercise_table) <- ExerciseIDtemp$V2

length(which(rownames(Exercise_table) == ExerciseIDtemp$V2))


#sanity checks
rownames(Exercise_ID) <- Exercise_ID$V1
rownames(ExerciseIDtemp) <- ExerciseIDtemp$V1
#sanity check to make sure that nothing was switched around.
Exercise_ID["0000704f7bd1dbca1898c70e12215e7e3a197745", "V2"] == ExerciseIDtemp["0000704f7bd1dbca1898c70e12215e7e3a197745", "V2"]

Exercise_ID["008ce2ba87ccc9fcf69117534438441c4afab4d8", "V2"] == ExerciseIDtemp["008ce2ba87ccc9fcf69117534438441c4afab4d8", "V2"]

#looks good.

Exercise_table_long <- t(Exercise_table)
dim(Exercise_table_long)
Exercise_taxa <- assignTaxonomy(Exercise_table_long, "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", multithread = 20)
dim(Exercise_taxa)

saveRDS(Exercise_taxa, file="tax_final.rds")
saveRDS(Exercise_table_long, file="seqtab_final.rds")


