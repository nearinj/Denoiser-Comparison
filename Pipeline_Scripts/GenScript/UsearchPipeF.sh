#!/bin/bash

#vars4 = sample list 
echo "************************************************** Running Unoise3 Pipeline *************************************************"

echo "Reading in " $1

readarray -t vars <$1

if [ ${vars[0]} -eq 1 ]; then
    echo "Removing old data"
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/Unoise
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/filtered_fastqs/rename
fi


cd ../../../DenoiseCompare_Out/${vars[1]}/$2/

mkdir Unoise

cd Unoise

echo "************************************ Pooling Samples into one file ************************************"
for Sample in ${vars[4]}
do
    echo $Sample
    ../../../../DenoiseCompare/Pipeline_Scripts/Medium_Stringent/rename.py -f ../filtered_fastqs/$Sample"_R1_filt.fastq" -s $Sample -o ../filtered_fastqs/rename
    cat ../filtered_fastqs/rename/$Sample"_R1_rename.fastq" >> all.merged.fq
done

	      
source activate qiime1


echo "********************************** Convert to FASTA ****************************************"
run_fastq_to_fasta.pl -p ${vars[3]} -o fasta_files all.merged.fq

echo "************************************* Finding Unique Reads and Abundances ***********************"
(time usearch10 -fastx_uniques fasta_files/*fasta -sizeout -relabel Uniq -fastaout uniques.fa) 2>Unitime.txt

echo "******************************** DeNoise Sequences ****************************************************************"
(time usearch10 -unoise3 uniques.fa -zotus zotus.fa) 2>Dtime.txt

cp zotus.fa zotusfix.fa

#change zotu label to otu ---> bug in usearch10
#some reason causes biom file to not be recognized by biom .... need to figure out whats going on here.
sed -i -e 's/Zotu/Otu/g' zotus.fa

echo "******************************************* Making Otu Tables *****************************************************"
usearch10 -otutab all.merged.fq -zotus zotus.fa -biomout otutab.biom -mapout map.txt -otutabout otutab.txt

echo "********************************************* Rarify OTU Table ******************************************************"
biom summarize-table -i otutab.biom -o seqtab_summary.txt

rare="$(cat seqtab_summary.txt | awk 'NR >= 16 { print }' | awk -F" " '{print $2}' | head -1)"
rarefix=${rare/.0/}
single_rarefaction.py -i otutab.biom -o seqtab_rare.biom -d $rarefix

