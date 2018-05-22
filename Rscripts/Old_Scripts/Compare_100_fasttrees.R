#Blueberry 100 trees analysis


#DadaVopenRef directory
setwd("projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/final_combined/sequences/pynast_aligned_seqs/sep_aligned/DadaVOpen-ref/shuffles/distances_trees/")


DadaVOpenRef_DMs <- list()

test <- list()
i <- 0
test[[1]] <- read.table(paste("shuffle",i,"/weighted_unifrac_rare_dadaVopen-ref.txt",sep=""), sep="\t", check.names = F, stringsAsFactors = F, row.names = 1, header=T)

while(i < 100 ){
  DadaVOpenRef_DMs[[i+1]] <- read.table(paste("shuffle",i,"/weighted_unifrac_rare_dadaVopen-ref.txt",sep=""), sep="\t", check.names = F, stringsAsFactors = F, row.names = 1, header=T)
  i <- i + 1
}


format_DMS <- function(DM, row_prefix, col_prefix){
  ret_tab <- DM[grepl(row_prefix, rownames(DM)), grepl(col_prefix, colnames(DM))]
  
  ret_tab <- ret_tab[gsub(col_prefix, row_prefix, colnames(ret_tab)),]
  
  if(identical(gsub(row_prefix, "", rownames(ret_tab)), gsub(col_prefix, "", colnames(ret_tab)))){
    return(ret_tab)
  }else
    return("Error")
}


#for each DM in the list set rows to Dada and cols to Open
DadaVOpenRef_format <- list()
j <- 1
for(i in DadaVOpenRef_DMs){
  DadaVOpenRef_format[[j]] <- as.matrix(format_DMS(i, "Dada_", "open-ref_"))
  j <- j+1
}

get_mean_diag <- function(list){
  #list that will store all the diagonals
  diag_list <- list()
  j <- 1
  while(j < 101 ){
    diag_list[[j]] <- diag(list[[j]])
    j <- j +1
  }
  #turn list into matrix and compute stats
  
  #each row represents a sample
  matrix <- do.call("cbind", diag_list)
  
  #generate a boxplot that 
  average_dists <- apply(matrix, 2, function(x) mean(x))
  
  return(average_dists)
}

DadaVOpenRef_mean_diag <- get_mean_diag(DadaVOpenRef_format)
hist(DadaVOpenRef_mean_diag)
boxplot(DadaVOpenRef_mean_diag)





#lets load this up for DadaVUnoise


setwd("../../../DadavUnoise/shuffles/distances_trees/")

DadaVUnoise_DMs <- list()


#read in DMs
j <- 0
while(j < 100){
  DadaVUnoise_DMs[[j+1]] <- read.table(paste("shuffle",j,"/weighted_unifrac_rare_DadavUnoise.txt",sep=""), sep="\t", check.names = F, stringsAsFactors = F, row.names = 1, header=T)
  j <- j +1
}


#format DMs
DadaVunoise_format <- list()
j <- 1
for(i in DadaVUnoise_DMs){
  DadaVunoise_format[[j]] <- as.matrix(format_DMS(i, "Dada_", "Unoise_"))
  j <- j+1
}

# get mean diags
DadaVUnoise_mean_diags <- get_mean_diag(DadaVunoise_format)
hist(DadaVUnoise_mean_diags)




#DeblurVDada

setwd("../../../DeblurVDada/shuffles/distance_trees/")

DeblurVDada_Dms <- list()

j <- 0
while(j < 100){
  DeblurVDada_Dms[[j+1]] <- read.table(paste("shuffle",j,"/weighted_unifrac_rare_DeblurVDada.txt",sep=""), sep="\t", check.names = F, stringsAsFactors = F, row.names = 1, header=T)
  j <- j +1
}

#format DMs
DeblurVDada_format <- list()
j <- 1
for(i in DeblurVDada_Dms){
  DeblurVDada_format[[j]] <- as.matrix(format_DMS(i, "Deblur_", "Dada_"))
  j <- j+1
}
DeblurVDada_avg_dist <- get_mean_diag(DeblurVDada_format)
hist(DeblurVDada_avg_dist)



#most likely not worth making 100 trees for each analysis as the differences are very small.
