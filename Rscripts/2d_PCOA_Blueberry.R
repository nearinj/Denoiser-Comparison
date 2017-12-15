#make 2d PCoA plots for blueberry


setwd("~/projects/DenoiseCompare_Out/Blueberry/med/COMBINED/")
#read in proportions for the PC
proportions <- read.table("plots/bdiv-out/weighted_unifrac_pc.txt", sep = "\t", nrow=1, skip=4)
#read in the number oof samples
SampleNum <- read.table("plots/bdiv-out/weighted_unifrac_pc.txt", sep = "\t", nrow=1)

#set cords for samplles

sample_cord <- read.table("plots/bdiv-out/weighted_unifrac_pc.txt", 
                               sep ="\t", skip = 9, 
                               nrow = SampleNum[[2]],
                               header=FALSE, row.names=1)
Shapes <-rep(x = 0, nrow(sample_cord))
Colors <- rep(x = "green", nrow(sample_cord))

Shapes[grep("Dada", rownames(sample_cord)) ] <- 15
Shapes[grep("Unoise", rownames(sample_cord))] <- 16
Shapes[grep("Deblur", rownames(sample_cord))] <- 17

Rhizosphere <- c("Bact3z2z2R", "Bact1z3z2R", "Bact198", "Bact183", "Bact181", "Bact197", 
                 "Bact5z3z2R", "Bact184", "Bact182", "Bact5z2z2R", "Bact196", "Bact1z1z2R",
                 "Bact1z4z2R", "Bact1z5x2R", "Bact3z4z2R", "Bact3z3z2R", "Bact5z4z2R",
                 "Bact3z5z2R", "Bact5z1z2R", "Bact180", "Bact195", "Bact1z2z2R", "Bact5z5z2R", "Bact3z1z2R")

for (val in Rhizosphere){
  Colors[grep(val, rownames(sample_cord))] <- "red"
}

paste(Colors)
plot(set1_sample_cord$V2, set1_sample_cord$V3, 
     xlab=paste("PC1 (",round(proportions$V1, 2)*100,"%)", sep=""),  
     ylab=paste("PC2 (",round(proportions$V2, 2)*100,"%)", sep=""),
     pch=c(Shapes), bg=Colors, col=Colors)
