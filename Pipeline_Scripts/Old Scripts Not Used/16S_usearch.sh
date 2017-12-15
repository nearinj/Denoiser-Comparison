#!/bin/bash

source activate dada2
echo "**************************************** FOR 16S DATA  ******************************************"

read -p "clear previous data?: " op

if [ $op -eq 1 ]; then
    echo "removing Old Information"
    rm -rf ../../../DenoiseCompare_Out/MediumUnoise3/
fi



echo "*************************** Running UNoise with Medium Strigencies *****************************************"

#read in the number of samples
read -p "Number of Paired Samples: " sNum


echo "**************************** Cutting Primers *********************************************************"

#create directories for the project
mkdir ../../../DenoiseCompare_Out/MediumUnoise3

mkdir ../../../DenoiseCompare_Out/MediumUnoise3/HMP


cd ../../../DenoiseCompare_Out/MediumUnoise3/HMP


mkdir primer_trimmed_fastqs


parallel --link --jobs $sNum --eta \
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



echo "******************************************* Filtering Reads ****************************************************************"
dada2_filter.R -f primer_trimmed_fastqs --truncLen 270,210 --maxN 0 --MaxEE 3,3 --truncQ 2 --threads $sNum --f_match _R1.*fastq.gz --r_match _R2.fastq.gz


echo "********************************************************* Pool Samples together ****************************************"
for Sample in Run1a Run2b Run3c Run4
do
    usearch10 -fastq_mergepairs ~/Mock_Communities/MockData/HMP/V4-V5-Rename/${Sample}*_R1.fastq -fastqout $Sample.merged.fq -relabel $Sample.
    cat $Sample.merged.fq >> all.merged.fq
done

source activate qiime1

echo "***************************************************** Convert to FASTA ***********************************************************"
run_fastq_to_fasta.pl -p $sNum -o fasta_files all.merged.fq


echo "**************************************************** Finding Unique Reads and Abundances *******************************************"
usearch10 -fastx_uniques fasta_files/*fasta -sizeout -relabel Uniq -fastaout uniques.fa

echo "******************************** DeNoise Sequences ****************************************************************"
usearch10 -unoise3 uniques.fa -zotus zotus.fa



#change zotu label to otu ---> bug in usearch10
#some reason causes biom file to not be recognized by biom .... need to figure out whats going on here.
sed -i -e 's/Zotu[0-9]/Otu/g' zotus.fa

echo "******************************************* Making Otu Tables *****************************************************"
usearch10 -otutab all.merged.fq -zotus zotus.fa -biomout otutab.biom -mapout map.txt -otutabout otutab.txt

echo "********************************************* Rarify OTU Table ******************************************************"
biom summarize-table -i otutab.biom -o seqtab_summary.txt

rare="$(cat seqtab_summary.txt | awk 'NR >= 16 { print }' | awk -F" " '{print $2}' | sort -n | head -1)"
rarefix=${rare/.0/}
single_rarefaction.py -i otutab.biom -o seqtab_rare.biom -d $rarefix





 
