### About
Folders contain the pipeline config files for the Trim.sh script that can be found in the Pipelines/Filter Scripts directory. The Old_Config folder contains configuration folders for an older data set that was not used in the final analysis for this study.


### About the Config Files
All study folders contain three different configuration files named High.txt, Med.txt, and Low.txt. There names are in accordance with the strigency to be used when filtering the fastq files. Each new line in the config file is taken in as a variable in the Trim.sh script the order that it is seen (i.e line 1 = variable $1). 

Line 1: Legacy variable that is no longer used
Line 2: When set to 1 will delete all previous data that was generated from running the Trim.sh script with that config file.
Line 3: The Directory name of where the raw fastq files are for the sample. Note (Trim.sh assumes all raw files are in ~/Mock_Communitites/MockData/ ). 
Line 4: The max amount of jobs/threads to use when running the Trim.sh script
Line 5: The filter stringency that the Trim.sh script should be run at.
Line 6: The forward primer sequence
Line 7: The Reverse primer sequence
Line 8: The length that the forward reads should be trimed to
Line 9: The length that the reverse reads should be trimed to
