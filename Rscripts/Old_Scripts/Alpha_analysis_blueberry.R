#alpha_analysis_blueberry

library(ggplot2)
library(reshape2)
library(ggthemes)
library(cowplot)

#functions

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
  samplenames <- gsub("Unoise_", "", samplenames)
  samplenames <- gsub("Deblur_", "", samplenames)
  samplenames <- gsub("open-ref_", "", samplenames)
  samplenames <- gsub("Deblur-noPos_", "", samplenames)
  tab$sample_names <- samplenames
  
  return(tab)
}


#load in the alpha data

setwd("~/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/biom/fixed_combined/")


alpha_data <- read.table("rare_alpha_div.tsv",  header=T, row.names=1, sep="\t", stringsAsFactors = F, check.names = F)

Pipes <- set_pipeline(alpha_data)

alpha_data$Pipelines <- Pipes

alpha_data <- add_sample_variable(alpha_data)


by_sample <- data.frame(Samples=character(64),
                        DADA2_Shannon=double(64),
                        Deblur_Shannon=double(64),
                        UNOISE3_Shannon=double(64),
                        DeblurnoPos_Shannon=double(64),
                        Openref_Shannon=double(64),
                        DADA2_Observed_OTU=integer(64),
                        Deblur_Observed_OTU=integer(64),
                        UNOISE3_Observed_OTU=integer(64),
                        DeblurnoPos_Observed_OTU=integer(64),
                        Openref_Observed_OTU=integer(64))


samples_to_compare <- alpha_data[alpha_data$Pipelines=="DADA2", "sample_names"]
alpha_data_fix <- alpha_data[alpha_data$sample_names %in% samples_to_compare,]
open_ref_keeped <- alpha_data[alpha_data_fix$Pipelines=="Open-Ref", "sample_names"]

samples_to_compare == open_ref_keeped
alpha_data_open <- alpha_data_fix[alpha_data_fix$Pipelines=="Open-Ref",]
alpha_data_open <- alpha_data_open[match(samples_to_compare, alpha_data_open$sample_names),]
open_ref_keeped <- alpha_data_open[,"sample_names"]

alpha_data_deblur <- alpha_data[alpha_data$Pipelines=="Deblur",]
alpha_data_deblur <- alpha_data_deblur[match(samples_to_compare, alpha_data_deblur$sample_names),]

samples_to_compare == open_ref_keeped

by_sample$Samples <- samples_to_compare
#samples_to_compare <- lapply(samples_to_compare, function(x) paste("open-ref_",x,sep=""))

by_sample$DADA2_Shannon <- alpha_data[alpha_data$Pipelines=="DADA2", "shannon"]
by_sample$Deblur_Shannon <- alpha_data_deblur[alpha_data_deblur$Pipelines=="Deblur", "shannon"]
by_sample$UNOISE3_Shannon <- alpha_data[alpha_data$Pipelines=="UNOISE3", "shannon"]
by_sample$Openref_Shannon <- alpha_data_open[alpha_data_open$Pipelines=="Open-Ref", "shannon"]
#by_sample$Openref_Shannon <- alpha_data[rownames(alpha_data) %in% samples_to_compare, "shannon"]

by_sample$DADA2_Observed_OTU <- alpha_data[alpha_data$Pipelines=="DADA2", "observed_otus"]
by_sample$Deblur_Observed_OTU <- alpha_data_deblur[alpha_data_deblur$Pipelines=="Deblur", "observed_otus"]
by_sample$UNOISE3_Observed_OTU <- alpha_data[alpha_data$Pipelines=="UNOISE3", "observed_otus"]
by_sample$Openref_Observed_OTU <- alpha_data_open[alpha_data_open$Pipelines=="Open-Ref", "observed_otus"]
#by_sample$Openref_Observed_OTU <- alpha_data[rownames(alpha_data) %in% samples_to_compare, "observed_otus"]
#sanity checks

#error checking to make sure by_sample is in the correct order.
j <- 0
for(i in samples_to_compare){

  if(!identical(alpha_data[paste("Dada_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"DADA2_Shannon"])){
    message("error")
    break
  }else
    j <- j +1
  
  if(!identical(alpha_data[paste("Deblur_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"Deblur_Shannon"])){
    message("error")
    break
  }else
    j <- j +1
  
  if(!identical(alpha_data[paste("Unoise_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"UNOISE3_Shannon"])){
    message("error")
    break
  }else
    j <- j +1
  
  if(!identical(alpha_data[paste("open-ref_",i,sep=""),"shannon"], by_sample[by_sample$Samples==i,"Openref_Shannon"])){
    message(alpha_data[paste("open-ref_",i,sep=""),"shannon"])
    message(by_sample[by_sample$Samples==i,"Openref_Shannon"])
    break
  }else
    j <- j + 1
  ##### check otu
  if(!identical(alpha_data[paste("Dada_",i,sep=""),"observed_otus"], by_sample[by_sample$Samples==i,"DADA2_Observed_OTU"])){
    message("error")
    break
  }else
    j <- j +1
  
  if(!identical(alpha_data[paste("Deblur_",i,sep=""),"observed_otus"], by_sample[by_sample$Samples==i,"Deblur_Observed_OTU"])){
    message("error")
    break
  }else
    j <- j +1
  
  if(!identical(alpha_data[paste("Unoise_",i,sep=""),"observed_otus"], by_sample[by_sample$Samples==i,"UNOISE3_Observed_OTU"])){
    message("error")
    break
  }else
    j <- j +1
  
  if(!identical(alpha_data[paste("open-ref_",i,sep=""),"observed_otus"], by_sample[by_sample$Samples==i,"Openref_Observed_OTU"])){
    message("error")
    break
  }else
    j <- j + 1
}
paste("number of passes",j,sep=" ")

#add linear regression line and equation to each graph
# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: http://goo.gl/K4yh

lm_eqn <- function(df, x ,y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

#plot 
test<-lm(DADA2_Shannon ~ Deblur_Shannon, by_sample)

DADAvDeblur_Shannon <- ggplot(by_sample, aes(x=DADA2_Shannon, y=Deblur_Shannon)) + geom_point() + theme_gdocs() +xlim(7.5,11) + ylim(7.5,11) +
  geom_smooth(method='lm')+ 
  annotate("text", x = 9.5, y = 10.75, label = lm_eqn(by_sample, by_sample$DADA2_Shannon, by_sample$Deblur_Shannon), colour="black", size = 5, parse=TRUE)

DADAvDeblur_Shannon

DadavUnoise_Shannon <- ggplot(by_sample, aes(x=DADA2_Shannon, y=UNOISE3_Shannon)) + geom_point() + theme_gdocs() + xlim(7.5,11) + ylim(7.5,11) +
  geom_smooth(method='lm')+
  annotate("text", x = 9.5, y = 10.75, label = lm_eqn(by_sample, by_sample$DADA2_Shannon, by_sample$UNOISE3_Shannon), colour="black", size = 5, parse=TRUE)
DadavUnoise_Shannon

UnoiseVDeblur_Shannon <- ggplot(by_sample, aes(x=UNOISE3_Shannon, y=Deblur_Shannon)) + geom_point() + theme_gdocs()  + xlim(7.5,11) + ylim(7.5,11) +
  geom_smooth(method='lm')+
  annotate("text", x = 9.5, y = 10.75, label = lm_eqn(by_sample, by_sample$UNOISE3_Shannon, by_sample$Deblur_Shannon), colour="black", size = 5, parse=TRUE)
UnoiseVDeblur_Shannon

DADAvDeblur_OTU <- ggplot(by_sample, aes(x=DADA2_Observed_OTU, y=Deblur_Observed_OTU)) + geom_point() + theme_gdocs() + xlim(400,2500) +ylim(400,2500) +
  geom_smooth(method='lm')+
  annotate("text", x = 1250, y = 2250, label = lm_eqn(by_sample, by_sample$DADA2_Observed_OTU, by_sample$Deblur_Observed_OTU), colour="black", size = 5, parse=TRUE)
DADAvDeblur_OTU

DadavUnoise_OTU <- ggplot(by_sample, aes(x=DADA2_Observed_OTU, y=UNOISE3_Observed_OTU)) + geom_point() + theme_gdocs() + xlim(400,2500) +ylim(400,2500) +
  geom_smooth(method='lm')+
  annotate("text", x = 1250, y = 2250, label = lm_eqn(by_sample, by_sample$DADA2_Observed_OTU, by_sample$UNOISE3_Observed_OTU), colour="black", size = 5, parse=TRUE)
DadavUnoise_OTU

UnoisevDeblur_OTU <- ggplot(by_sample, aes(x=UNOISE3_Observed_OTU, y=Deblur_Observed_OTU)) + geom_point() + theme_gdocs() + xlim(400,2500) +ylim(400,2500) +
  geom_smooth(method='lm')+
  annotate("text", x = 1250, y = 2250, label = lm_eqn(by_sample, by_sample$UNOISE3_Observed_OTU, by_sample$Deblur_Observed_OTU), colour="black", size = 5, parse=TRUE)
UnoisevDeblur_OTU

DeblurVOpen_shannon <- ggplot(by_sample, aes(x=Deblur_Shannon, y=Openref_Shannon)) +geom_point() + theme_gdocs()  + xlim(7.5,11) + ylim(7.5,11) +
  geom_smooth(method='lm')+
  annotate("text", x = 9.5, y = 10.75, label = lm_eqn(by_sample, by_sample$Deblur_Shannon, by_sample$Openref_Shannon), colour="black", size = 5, parse=TRUE)
DeblurVOpen_shannon

UnoiseVOpen_shannon <- ggplot(by_sample, aes(x=UNOISE3_Shannon, y=Openref_Shannon)) + geom_point() + theme_gdocs()  + xlim(7.5,11) + ylim(7.5,11)+
  geom_smooth(method='lm')+
  annotate("text", x = 9.5, y = 10.75, label = lm_eqn(by_sample, by_sample$UNOISE3_Shannon, by_sample$Openref_Shannon), colour="black", size = 5, parse=TRUE)
UnoiseVOpen_shannon

DadaVOpen_shannon <- ggplot(by_sample, aes(x=DADA2_Shannon, y=Openref_Shannon)) + geom_point() + theme_gdocs()  + xlim(7.5,11) + ylim(7.5,11)+
  geom_smooth(method='lm')+
  annotate("text", x = 9.5, y = 10.75, label = lm_eqn(by_sample, by_sample$DADA2_Shannon, by_sample$Openref_Shannon), colour="black", size = 5, parse=TRUE)
DadaVOpen_shannon

DeblurVOpen_OTU <- ggplot(by_sample, aes(x=Deblur_Observed_OTU, y=Openref_Observed_OTU)) + geom_point() + theme_gdocs()  + xlim(400,2500) +ylim(400,2500)+
  geom_smooth(method='lm')+
  annotate("text", x = 1250, y = 2250, label = lm_eqn(by_sample, by_sample$Deblur_Observed_OTU, by_sample$Openref_Observed_OTU), colour="black", size = 5, parse=TRUE)
DeblurVOpen_OTU

UnoiseVOpen_OTU <- ggplot(by_sample, aes(x=UNOISE3_Observed_OTU, y=Openref_Observed_OTU)) + geom_point() + theme_gdocs()  + xlim(400,2500) +ylim(400,2500)+
  geom_smooth(method='lm')+
  annotate("text", x = 1250, y = 2250, label = lm_eqn(by_sample, by_sample$UNOISE3_Observed_OTU, by_sample$Openref_Observed_OTU), colour="black", size = 5, parse=TRUE)
UnoiseVOpen_OTU

DadaVOpen_OTU <- ggplot(by_sample, aes(x=DADA2_Observed_OTU, y=Openref_Observed_OTU)) + geom_point() + theme_gdocs() + xlim(400,2500) +ylim(400,2500)+
  geom_smooth(method='lm')+
  annotate("text", x = 1250, y = 2250, label = lm_eqn(by_sample, by_sample$DADA2_Observed_OTU, by_sample$Openref_Observed_OTU), colour="black", size = 5, parse=TRUE)
DadaVOpen_OTU

#set working directory to new save location

setwd("~/projects/DenoiseCompare/Revised_figures/")
#Alpha Shannon plot


Shannon_grid <- plot_grid(DadavUnoise_Shannon, DADAvDeblur_Shannon, UnoiseVDeblur_Shannon, DadaVOpen_shannon, DeblurVOpen_shannon, UnoiseVOpen_shannon, labels="AUTO")
Shannon_grid

#Alpha OTU plot

OTU_grid <- plot_grid(DadavUnoise_OTU, DADAvDeblur_OTU, UnoisevDeblur_OTU, DadaVOpen_OTU, DeblurVOpen_OTU, UnoiseVOpen_OTU, labels="AUTO")
OTU_grid
