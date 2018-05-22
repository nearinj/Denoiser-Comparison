#!/bin/bash

#vars4 = sample list 
echo "************************************************** Running Unoise3 Pipeline *************************************************"

echo "Reading in " $1

readarray -t vars <$1

if [ ${vars[0]} -eq 1 ]; then
    echo "Removing old data"
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/Unoise
fi


cd ../../../DenoiseCompare_Out/${vars[1]}/$2/

if [[ ${vars[6]} != "" ]]; then
    echo "running on "${vars[6]}

       cd ${vars[6]} 
       mkdir Unoise
       cd Unoise
       for i in ../*R1*; do
	   name=${i/_*/}
	   usearch10 -fastq_mergepairs ${name}*_R1*.fastq -fastqout $name.merged.fq -relabel $name
	   cat $name.merged.fq >> all.merged.fq
       done
elif [[ ${vars[7]} != "" ]]; then
	 echo "running on" ${vars[7]}

	 mkdir Unoise

	 cd Unoise
	 mkdir Out
	 for i in ../filtered_fastqs/*R1*; do
	     prefix="../filtered_fastqs/"
	     name=${i#$prefix}
	     usearch10 -fastq_mergepairs $i -fastqout Out/$name.merged.fq -relabel $name
	     cat Out/$name.merged.fq >> all.merged.fq
	 done
	 
    
else   
    mkdir Unoise

    cd Unoise

    pwd

    echo "***********************************Pooling Samples Together ***********************************"
    for Sample in ${vars[4]}
    do
	usearch10 -fastq_mergepairs ../filtered_fastqs/${Sample}*_R1*.fastq -fastqout $Sample.merged.fq -relabel $Sample.
	cat $Sample.merged.fq >> all.merged.fq
    done
fi

source activate qiime1

echo "********************************** Convert to FASTA ****************************************"
run_fastq_to_fasta.pl -p ${vars[3]} -o fasta_files all.merged.fq

echo "************************************* Finding Unique Reads and Abundances ***********************"
(/usr/bin/time -v usearch10 -fastx_uniques fasta_files/*fasta -sizeout -relabel Uniq -fastaout uniques.fa) 2>Unitime.txt

echo "******************************** DeNoise Sequences ****************************************************************"
(/usr/bin/time -v usearch10 -unoise3 uniques.fa -zotus zotus.fa) 2>Dtime.txt

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
