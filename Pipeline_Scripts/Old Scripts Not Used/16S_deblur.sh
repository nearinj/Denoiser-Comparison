#!/bin/bash


echo "************************************************** FOR 16S DATA  ********************************************************"
source activate dada2
read -p "clear previous data? " op

if [ $op -eq 1 ]; then
    echo "removing"
    rm -rf ../../../DenoiseCompare_Out/MediumDeblur/HMP
fi

echo "******************************************************************* RUNING DEBLUR WITH MEDIUM STRIGENCIES ***********************"

read -p "Number of Paired Samples: " sNum


echo "*********************** Cutting Primers **********************************************"

mkdir ../../../DenoiseCompare_Out/MediumDeblur

mkdir ../../../DenoiseCompare_Out/MediumDeblur/HMP

cd ../../../DenoiseCompare_Out/MediumDeblur/HMP


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

#Filter the reads
dada2_filter.R -f primer_trimmed_fastqs --truncLen 270,210 --maxN 0 --maxEE 3,3 --truncQ 2 --threads $sNum --f_match _R1.*fastq.gz --r_match _R2.fastq.gz

echo "******************************** pair samples ***********************************"
source activate qiime1
run_pear.pl -p $sNum -o stiched_reads primer_trimmed_fastqs/*fastq.gz

mkdir paired

mv stiched_reads/*.assembled.fastq paired



source activate deblurenv
echo "**************************************************** Run Through Deblur workflow ***********************************************"
deblur workflow --seqs-fp paired --output-dir output -t 250 -a 10

echo "******************************************************Summarize biom file and Rarifying to highest depth possible **********************"
source activate qiime1

biom summarize-table -i output/all.biom -o summary.txt

#get lowest number of reads in a sample and rarify
rare="$(cat summary.txt | awk 'NR >= 16 { print }' | awk -F" " '{print $2}' |head -1)"
rarefix=${rare/.0/}
single_rarefaction.py -i output/all.biom -o output/all.rare.biom -d $rarefix
