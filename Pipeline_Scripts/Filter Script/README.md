### Rarefy_FiltReads.sh

Script used to generate the rarefied reads used to measure the speed and memory usage of the three different pipelines. This script takes the same config file that was used to generate the filtered reads using the TrimFilt.sh script.

### TrimFilt.sh

Script that trims off primers from raw fastq files and quaility filters them depending on the config file that is provided. The config files can be found in the Config/Trim directory.

### demult_barcode_readnames.py

Script that takes in a raw gzipped fastq file, a metadata file containing the sample names and the corresponding barcodes and demultiplexes the reads into individual fastq files for each sample found in the metadata file.

### rename_fastq_samples.py
Script that takes in a fastq file and a sample name and then renames the fastq samples so that they corporate with all three pipelines.
