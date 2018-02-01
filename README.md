## About
The aim of this project is to compare the results from three different Illumina sequence denoising pipelines.
  
The three pipelines are:
* [DADA2](https://benjjneb.github.io/dada2/)
* [UNOISE3](https://www.drive5.com/usearch/manual/cmd_unoise3.html)
* [Deblur](https://github.com/biocore/deblur)

None of the raw sequencing data is contained within this repository.

### Blast_db
Contains all scripts used to modify any of the Database files
All Databases were generated using blast tools

### Config
Cotains all the configuration files to run the different pipelines

### Data_Analysis
Contains all scripts that were used strictly for the analysis of the resulting amplicon sequence variants from the pipelines

### Pipeline_Scripts
Contains the scripts that were used to run the three different pipelines

### Rscripts
Contains all the Rscripts that were used to generate plots as well as assign taxonomy.
