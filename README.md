# Denoiser-Comparison
## About
   The aim of this project is to compare the results from three different Illumina 454 sequence denoising pipelines.
   The three pipelines are:
       DADA2
       UNOISE3
       Deblur
   None of the raw data is contained within this repository as it would easily exceed the maximum limit
### Blast_db
Contains all scripts used to modify any of the Database files
All Databases were generated using blast tools
slice_amplified_region.py - Takes a fasta file and two primers and slices out the region found between the two primers in the fasta sequences
slice_amplified_region_forward-only_set-length.py - Takes a fasta file and one forward primer as well as a base pair length. It slices out the region from the forward parimer to the end of the base pair length.
### Config
Cotains all the configuration files to run the different pipelines
### Data_Analysis
Contains all scripts that were used stricktly for the analysis of the resulting ASV from the pipelines
### Pipeline_Scripts
Contains the scripts that were used to run the three different pipelines
### Rscripts
Contains all the Rscripts that were used to generate plots as well as assign taxonomy.
