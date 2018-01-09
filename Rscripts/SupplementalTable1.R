#Run after New_Abudance_Stacked_Bars.R

#Used to create Supplmentary table 1

setwd("~/projects/DenoiseCompare/Final_Figures/")
HMP_Table <- HMP_MED[-c(1,20), -c(2,3)]

HMP_Table <- findPresence(HMP_Table)
write.table(HMP_Table, file="HMP_Organism_Presence.tsv", quote = F, sep="\t", col.names= NA)

Zymock_Table <- Zymock_MED[-c(1,10), -c(2,3)]
Zymock_Table <- findPresence(Zymock_Table)
write.table(Zymock_Table, file="Zymock_Organism_Presence.tsv", quote=F, sep="\t", col.names=NA)

Extreme_Table <- Mock12_MED[-c(1,29), -c(2,3)]
Extreme_Table <- findPresence(Extreme_Table)
write.table(Extreme_Table, file="Extreme_Organism_Presence.tsv", quote=F, sep="\t", col.names=NA)

Fungal_Table <- Mock9_MED[-c(1,15), -c(2,3)]
Fungal_Table <- findPresence(Fungal_Table)
write.table(Fungal_Table, file="Fungal_Organism_Presence.tsv", quote=F, sep="\t", col.names=NA)

findPresence <- function(tab) {
  tab$Dada[tab$Dada > 0] <- "Yes"
  tab$Deblur[tab$Deblur > 0] <- "Yes"
  tab$Unoise[tab$Unoise > 0] <- "Yes"
  
  tab$Dada[tab$Dada == 0] <- "No"
  tab$Deblur[tab$Deblur == 0] <- "No"
  tab$Unoise[tab$Unoise == 0] <- "No"
  
  return(tab)
}
