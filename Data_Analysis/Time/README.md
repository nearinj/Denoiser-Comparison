### BlueBerry_Memory_Plot.R
Script used to generate the memory usage plot while running the Blueberry data at different amounts of rarefied reads for each pipeline.
### Blueberry_time_data.tsv
Table generated from running pull_time.py on each of the rarefied read directories (i.e Filt5k, Filt10k, ..... Filt40k).
### Time_Data_for_Blueberry.R
Script used to generate the time usage plot while running the Blueberry data at different amounts of rarefied reads for each pipeline.
### pull_time.py
Script that takes in a directory were there are different amounts of rarefied reads in multiple directories and extracts the user time as well as the peak amount of memory used by each pipeline while running those reads.
