### About
This folder contains all config files used to run the Dada2Pipe.sh script. Each File is named by the sample set that the config file corresponds with.

### Config File Setup

Each new line in the config file is taken in as a variable in the Dada2Pipe.sh script in the order that it is seen (i.e line 1 = variable $0).

Line 1: When set to 1 will delete all previous data generated by the Dada2Pipe.sh script that was run using this config file.

Line 2: The directory that contains the folder that contains the filtered reads.

Line 3: Legacy variable that is no longer used

Line 4: The max amount of threads to use

Line 6: When set to Single will run the DaDa2 inference script using unpaired mode

Line 7: When set to a specific directory will run the Dada2Pipeline on the reads within that directory (Used for processes the reads used to determine the speed and memory usage of the pipelines. See Blueberry_time directory)
