## About
The aim of this project is to compare the results from three different Illumina sequence denoising pipelines.
  
The three pipelines are:
* [DADA2](https://benjjneb.github.io/dada2/)
* [UNOISE3](https://www.drive5.com/usearch/manual/cmd_unoise3.html)
* [Deblur](https://github.com/biocore/deblur)

We also compared the pipelines to Open Reference OTU clustering at 97%

None of the raw sequencing data is contained within this repository.

### Blast_db
Contains all scripts used to modify any of the Database files
  More specifically it contains the script that was used to determine the number of unique sequences in each expected sequence file. 
All Databases were generated using blast tools

### Config
Cotains all the configuration files to run the different pipelines on each real dataset and mock dataset

### Data_Analysis
Contains all scripts that were used strictly for the analysis of the resulting amplicon sequence variants from the pipelines including both ASV type analysis and ASV abundance analysis.

### Pipeline_Scripts
Contains the scripts that were used to run the four different pipelines (DADA2, Deblur, UNOISE3 and a VSEARCH based open reference OTU clustering pipeline).

### Rscripts
Contains all the Rscripts that were used to generate plots for the manuscript as well as scripts that were used to assign taxonomy to biom tables from real datasets (i.e not mock communities). 
