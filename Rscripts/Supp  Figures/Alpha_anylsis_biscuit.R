#alpha diversity analysis for biscuit data...

library(ggplot2)
library(reshape2)
library(ggthemes)
library(cowplot)

setwd("~/projects/DenoiseCompare_Out/biscuit/med/COMBINED/biom/fixed_combined/")



alpha_data <- read.table("rare_alpha_div.tsv", header=T, row.names=1, sep="\t", stringsAsFactors = F, check.names = F)

set_pipeline <- function(tab){
  test <- rownames(alpha_data)
  pipes <- vector()
  for (i in 1:length(test)){
    if(grepl("Dada", test[i])){
      pipes[i] <- "DADA2"
    }else if(grepl("Deblur_",test[i])){
      pipes[i] <- "Deblur"
    }else if(grepl("Unoise",test[i])){
      pipes[i] <- "UNOISE3"
    }else if(grepl("Deblur-noPos_", test[i])){
      pipes[i] <- "Deblur-noPos"
    }else if(grepl("open-ref_", test[i])){
      pipes[i] <- "Open-Ref"
    }else
      pipes[i] <- "Other"
  }
  return(pipes)
}

#add sample names
add_sample_variable <- function(tab){
  samplenames <- rownames(tab)
  samplenames <- gsub("Dada_", "", samplenames)
  samplenames <- gsub("Unoise__", "", samplenames)
  samplenames <- gsub("Deblur_", "", samplenames)
  samplenames <- gsub("open-ref_", "", samplenames)
  samplenames <- gsub("Deblur-noPos_", "", samplenames)
  tab$sample_names <- samplenames
  
  return(tab)
}

Pipes <- set_pipeline(alpha_data)

#sanity check
length(which(Pipes == "DADA2")) == length(which(Pipes == "Deblur"))
length(which(Pipes == "UNOISE3")) == length(which(Pipes == "DADA2"))


alpha_data$Pipelines <- Pipes

alpha_data <- add_sample_variable(alpha_data)
alpha_data$sample_names <- gsub("001101", "", alpha_data$sample_names)



#clean up and graph number of observed OTUs against each other.

#graph DADA2 vs UNOISE3

by_sample <- data.frame(Samples=character(32),
                        DADA2_Shannon=double(32),
                        Deblur_Shannon=double(32),
                        UNOISE3_Shannon=double(32),
                        DeblurnoPos_Shannon=double(32),
                        Open_Shannon=double(32),
                        DADA2_Observed_ASV=integer(32),
                        Deblur_Observed_ASV=integer(32),
                        UNOISE3_Observed_ASV=integer(32),
                        DeblurnoPos_Observed_OTU=integer(32),
                        Open_Observed_OTU=integer(32))

samples_to_compare <- alpha_data[alpha_data$Pipelines=="DADA2", "sample_names"]
alpha_data_fix <- alpha_data[alpha_data$sample_names %in% samples_to_compare,]
open_ref_keeped <- alpha_data[alpha_data_fix$Pipelines=="Open-Ref", "sample_names"]

samples_to_compare == open_ref_keeped
alpha_data_open <- alpha_data_fix[alpha_data_fix$Pipelines=="Open-Ref",]
alpha_data_open <- alpha_data_open[match(samples_to_compare, alpha_data_open$sample_names),]
open_ref_keeped <- alpha_data_open[,"sample_names"]

samples_to_compare == open_ref_keeped

by_sample$Samples <- samples_to_compare
#samples_to_compare <- lapply(samples_to_compare, function(x) paste("open-ref_",x,sep=""))

by_sample$DADA2_Shannon <- alpha_data[alpha_data$Pipelines=="DADA2", "shannon"]
by_sample$Deblur_Shannon <- alpha_data[alpha_data$Pipelines=="Deblur", "shannon"]
by_sample$UNOISE3_Shannon <- alpha_data[alpha_data$Pipelines=="UNOISE3", "shannon"]
by_sample$Open_Shannon <- alpha_data_open[alpha_data_open$Pipelines=="Open-Ref", "shannon"]
#by_sample$Openref_Shannon <- alpha_data[rownames(alpha_data) %in% samples_to_compare, "shannon"]

by_sample$DADA2_Observed_ASV <- alpha_data[alpha_data$Pipelines=="DADA2", "observed_otus"]
by_sample$Deblur_Observed_ASV <- alpha_data[alpha_data$Pipelines=="Deblur", "observed_otus"]
by_sample$UNOISE3_Observed_ASV <- alpha_data[alpha_data$Pipelines=="UNOISE3", "observed_otus"]
by_sample$Open_Observed_OTU <- alpha_data_open[alpha_data_open$Pipelines=="Open-Ref", "observed_otus"]
#by_sample$Openref_Observed_OTU <- alpha_data[rownames(alpha_data) %in% samples_to_compare, "observed_otus"]


#error checking
for(i in samples_to_compare){
  
  if(!identical(alpha_data[paste("Dada_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"DADA2_Shannon"])){
    message("Dada doesn't match")
    break
  }
  
  if(!identical(alpha_data[paste("Deblur_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"Deblur_Shannon"])){
    message("Deblur doesn't match")
    break
  }
  
  if(!identical(alpha_data[paste("Unoise__",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"UNOISE3_Shannon"])){
    message("Unoise doesn't match")
    break
  }
  
  if(!identical(alpha_data[paste("open-ref_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"Open_Shannon"])){
    message("open-ref doesn't match")
    break
  }
  
}


#plot 
lm_eqn <- function(df, x ,y){
  m <- cor.test(y, x, method="spearman", exact=F);
  eq <- substitute(italic(r)~"="~r2*",  "~~italic(p)~"="~pv, list(r2 = format(m$estimate, digits = 3), pv = format(m$p.value, digits=3)))
  as.character(as.expression(eq));                 
}

#doesn't do the best at looking pairwise samples
DADAvDeblur_Shannon <- ggplot(by_sample, aes(x=DADA2_Shannon, y=Deblur_Shannon)) + geom_point() + xlim(2,7) + ylim(2,7) + theme_gdocs()+
  annotate("text", x = 3.5, y = 6.75, label = lm_eqn(by_sample, by_sample$DADA2_Shannon, by_sample$Deblur_Shannon), colour="black", size = 5, parse=TRUE)
DADAvDeblur_Shannon

DadavUnoise_Shannon <- ggplot(by_sample, aes(x=DADA2_Shannon, y=UNOISE3_Shannon)) + geom_point() + xlim(2,7) + ylim(2,7) + theme_gdocs()+
  annotate("text", x = 3.5, y = 6.75, label = lm_eqn(by_sample, by_sample$DADA2_Shannon, by_sample$UNOISE3_Shannon), colour="black", size = 5, parse=TRUE)
DadavUnoise_Shannon

UnoiseVDeblur_Shannon <- ggplot(by_sample, aes(x=UNOISE3_Shannon, y=Deblur_Shannon)) + geom_point() + xlim(2,7) + ylim(2,7) + theme_gdocs()+
  annotate("text", x = 3.5, y = 6.75, label = lm_eqn(by_sample, by_sample$UNOISE3_Shannon, by_sample$Deblur_Shannon), colour="black", size = 5, parse=TRUE)
UnoiseVDeblur_Shannon

DADAvDeblur_OTU <- ggplot(by_sample, aes(x=DADA2_Observed_ASV, y=Deblur_Observed_ASV)) + geom_point() + xlim(55,370) + ylim(55,370) + theme_gdocs()+
  annotate("text", x = 150, y = 350, label = lm_eqn(by_sample, by_sample$DADA2_Observed_ASV, by_sample$Deblur_Observed_ASV), colour="black", size = 5, parse=TRUE)
DADAvDeblur_OTU

DadavUnoise_OTU <- ggplot(by_sample, aes(x=DADA2_Observed_ASV, y=UNOISE3_Observed_ASV)) + geom_point() + xlim(55,370) + ylim(55,370) + theme_gdocs()+
  annotate("text", x = 150, y = 350, label = lm_eqn(by_sample, by_sample$DADA2_Observed_ASV, by_sample$UNOISE3_Observed_ASV), colour="black", size = 5, parse=TRUE)
DadavUnoise_OTU

UnoisevDeblur_OTU <- ggplot(by_sample, aes(x=UNOISE3_Observed_ASV, y=Deblur_Observed_ASV)) + geom_point() + xlim(55,370) + ylim(55,370) + theme_gdocs()+
  annotate("text", x = 150, y = 350, label = lm_eqn(by_sample, by_sample$UNOISE3_Observed_ASV, by_sample$Deblur_Observed_ASV), colour="black", size = 5, parse=TRUE)
UnoisevDeblur_OTU

DeblurVOpen_shannon <- ggplot(by_sample, aes(x=Deblur_Shannon, y=Open_Shannon)) +geom_point() + theme_gdocs() +xlim(2,7) + ylim(2,7) +
  annotate("text", x = 3.5, y = 6.75, label = lm_eqn(by_sample, by_sample$Deblur_Shannon, by_sample$Open_Shannon), colour="black", size = 5, parse=TRUE)
DeblurVOpen_shannon

UnoiseVOpen_shannon <- ggplot(by_sample, aes(x=UNOISE3_Shannon, y=Open_Shannon)) + geom_point() + theme_gdocs() +xlim(2,7) + ylim(2,7) +
  annotate("text", x = 3.5, y = 6.75, label = lm_eqn(by_sample, by_sample$UNOISE3_Shannon, by_sample$Open_Shannon), colour="black", size = 5, parse=TRUE)
UnoiseVOpen_shannon

DadaVOpen_shannon <- ggplot(by_sample, aes(x=DADA2_Shannon, y=Open_Shannon)) + geom_point() + theme_gdocs() +xlim(2,7) + ylim(2,7) + 
  annotate("text", x = 3.5, y = 6.75, label = lm_eqn(by_sample, by_sample$DADA2_Shannon, by_sample$Open_Shannon), colour="black", size = 5, parse=TRUE)
DadaVOpen_shannon

DeblurVOpen_OTU <- ggplot(by_sample, aes(x=Deblur_Observed_ASV, y=Open_Observed_OTU)) + geom_point() + theme_gdocs() + xlim(55,370) + ylim(55,370)+
  annotate("text", x = 150, y = 350, label = lm_eqn(by_sample, by_sample$Deblur_Observed_ASV, by_sample$Open_Observed_OTU), colour="black", size = 5, parse=TRUE)
DeblurVOpen_OTU

UnoiseVOpen_OTU <- ggplot(by_sample, aes(x=UNOISE3_Observed_ASV, y=Open_Observed_OTU)) + geom_point() + theme_gdocs()  + xlim(55,370) + ylim(55,370)+
  annotate("text", x = 150, y = 350, label = lm_eqn(by_sample, by_sample$UNOISE3_Observed_ASV, by_sample$Open_Observed_OTU), colour="black", size = 5, parse=TRUE)
UnoiseVOpen_OTU

DadaVOpen_OTU <- ggplot(by_sample, aes(x=DADA2_Observed_ASV, y=Open_Observed_OTU)) + geom_point() + theme_gdocs()  + xlim(55,370) + ylim(55,370)+
  annotate("text", x = 150, y = 350, label = lm_eqn(by_sample, by_sample$DADA2_Observed_ASV, by_sample$Open_Observed_OTU), colour="black", size = 5, parse=TRUE)
DadaVOpen_OTU

#Alpha Shannon plot

#set working directory to image saving dir
setwd("~/projects/DenoiseCompare/Revised_figures/")
Shannon_grid <- plot_grid(DadavUnoise_Shannon, DADAvDeblur_Shannon, UnoiseVDeblur_Shannon, DadaVOpen_shannon, DeblurVOpen_shannon, UnoiseVOpen_shannon, labels="AUTO")
Shannon_grid

#Alpha OTU plot

OTU_grid <- plot_grid(DadavUnoise_OTU, DADAvDeblur_OTU, UnoisevDeblur_OTU, DadaVOpen_OTU, DeblurVOpen_OTU, UnoiseVOpen_OTU, labels="AUTO")
OTU_grid

#interesting that unoise finds higher OTUs per sample but doesn't find the highest amount of OTUs in total...