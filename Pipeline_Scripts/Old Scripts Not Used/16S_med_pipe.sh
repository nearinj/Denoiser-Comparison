#!/bin/bash

echo "*******************FOR 16S DATA************************"

read -p "clear previous data? " op
echo "*******************************************Clear Previous Data*****************************************"

if [ $op -eq 1 ]; then
    echo "removing"
    rm -rf ../../../DenoiseCompare_Out/MediumDada2/
fi


echo "**************** Running Dada2 With Medium Strigencies ******************************"
read -p "Directory of 16S Data: " DataDir
read -p "Number of Paired Samples: " sNum


echo "********* Cutting Primers (Primers are cut for all Pipelines)***********"

#create directories for the project
mkdir ../../../DenoiseCompare_Out/MediumDada2

mkdir ../../../DenoiseCompare_Out/MediumDada2/HMP

cd ../../../DenoiseCompare_Out/MediumDada2/HMP

mkdir primer_trimmed_fastqs

parallel --link --jobs $sNum --eta\
	   'cutadapt \
    --pair-filter any \
    --no-indels \
    --discard-untrimmed \
    -g GTGYCAGCMGCCGCGGTAA \
    -G CCGYCAATTYMTTTRAGTTT \
    -o primer_trimmed_fastqs/{1/}.gz \
    -p primer_trimmed_fastqs/{2/}.gz \
    {1} {2} \
    > primer_trimmed_fastqs/{1/}_cutadapt_log.txt' \
	     ::: ~/Mock_Communities/MockData/HMP/V4-V5-Rename/*_R1*.fastq ::: ~/Mock_Communities/MockData/HMP/V4-V5-Rename/*_R2*.fastq



#put cut logs together
parse_cutadapt_logs.py -i primer_trimmed_fastqs/*log.txt





echo "*********Filtering Reads******************" 
dada2_filter.R -f primer_trimmed_fastqs --truncLen 270,210 --maxN 0 --maxEE 3,3 --truncQ 2 --threads $sNum --f_match _R1.*fastq.gz --r_match _R2.fastq.gz


echo "*********************Getting Error Rates*********************"

dada2_inference.R -f filtered_fastqs/ --seed 1995 -t $sNum --verbose


echo "*******************************Chimera Checking Assigning Taxonomy******************************"

dada2_chimera_taxa.R -i seqtab.rds -r ~/etc/databases/rdp_train_set_16.fa --skip_species -t $sNum

echo "**********************************Merging Log Files*********************************************"

merge_logfiles.R -i cutadapt_log.txt,dada2_filter_read_counts.txt,dada2_inferred_read_counts.txt,dada2_nonchimera_counts.txt \
		         -n cutadapt,filter,infer,chimera -o combined_log.txt

echo "*******************************Converting to BIOM and FASTA**************************************"
convert_dada2_out.R -i \
		    seqtab_final.rds -b seqtab.biom -f seqtab.fasta --taxa_in tax_final.rds


echo "************************************Summarize biom file and Rarifying to highest depth possible****************************"
biom summarize-table -i seqtab.biom -o seqtab_summary.txt


source activate qiime1

#get lowest number of reads in a sample and rarify
rare="$(cat seqtab_summary.txt | awk 'NR >= 16 { print }' | awk -F" " '{print $2}' | sort -n | head -1)"
rarefix=${rare/.0/}
single_rarefaction.py -i seqtab.biom -o seqtab_rare.biom -d $rarefix



