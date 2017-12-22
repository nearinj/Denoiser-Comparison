### Old_Scripts
Contains all of the outdated scripts that were either updated and can be found here or were found to provide no information
### Assign_taxa_deblur_unoise.R
Script that is takes the a tsv convert biom file from the Deblur and UNOISE3 pipelines along with the corresponding ASV files (fasta file) and assigns taxonomy to ASV.
### Blast_Results_With_Lines.R
Takes the results from blasting against the expected sequences of each mock community, followed by blasting any none hits against a database and plots them into a stacked barchart. It also takes the total number of possible unique regions that could have been sequenced from the mock community during sequencing and represents it as a dotted black line over each plot.
### Bray_curtis_blueberry.R
Takes the results of taxonomy assignment and generation of a Bray-Curtis distance matrix and plots the inter-sample distances between the three different pipelines in box plots. It also plots the intra-sample distance between two distinct sample groups (bulk and rhizosphere) as a base line for the variability between samples.
### New_Abundance_Stacked_Bars.R
Takes the results from all the 97% or greater blast hits from the expected sequence database for each mock community and plots out their abundances. 
### Weighted_Unifrac_Blueberry.R
Takes the distance matrix generated from a phylogenetic tree that was generated from a biom table complied of all three pipelines where sequences arrived by other pipelines were filled in as 0 in abundance. It then plots the inter-sample unifrac distances between each pipeline and displays them as box plots.  It also plots the intra-sample distance between two distinct sample groups (bulk and rhizosphere) as a base line for the variability between samples. 
### Blueberry_ordination.R
Takes the a Bray-Curtis distance matrix and creates an NMDS plot where pipelines are represented by different shapes and the same samples are connected by lines as well as color.
