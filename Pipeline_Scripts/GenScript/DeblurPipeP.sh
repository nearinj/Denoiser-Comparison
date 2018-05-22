#!/bin/bash

source activate qiime2-2017.9

echo "******************************** Reading in Config File *******************************************************"
echo "Reading in " $1

readarray -t vars <$1

if [ ${vars[0]} -eq 1 ]; then
    echo "Removing old data"
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/deblur
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/filtered_fastqs/renamedpaired
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/${vars[6]}/stitched_reads
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/${vars[6]}/deblurP
fi


cd ../../../DenoiseCompare_Out/${vars[1]}/$2



if [[ ${vars[6]} != "" ]]; then
    cd ${vars[6]}
    run_pear.pl -o stitched_reads ./* -p ${vars[4]}

    cd stitched_reads
    mkdir renamedpaired

    for i in *assembled.fastq; do name=${i/.*/}; cp $PWD/$i $PWD/renamedpaired/$name"_X_L001_R1_001.fastq"; done

    cd renamedpaired

    gzip *

    cd ../../


else

	echo "******************************* Pairing reads together **************************************************"
	run_pear.pl -o stitched_reads filtered_fastqs/*fastq -p ${vars[4]}
	echo "******************************* Creating links to import into qiime2 artifact ********************************************"

	cd stitched_reads

	mkdir renamedpaired

	for i in *assembled.fastq; do name=${i/.*/}; cp  $PWD/$i $PWD/renamedpaired/$name"_X_L001_R1_001.fastq"; done

	cd renamedpaired

	gzip *


	cd ../../

	echo "***************************** Importing filtered files into qiime2 artifact **********************************************"
fi
mkdir deblurP

qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path stitched_reads/renamedpaired --source-format CasavaOneEightSingleLanePerSampleDirFmt --output-path deblurP/filtered.qza

cd deblurP

echo "******************** Running Deblur on forward reads ********************************************"

(/usr/bin/time -v qiime deblur denoise-16S --p-jobs-to-start ${vars[4]} --i-demultiplexed-seqs filtered.qza --p-trim-length 350 --o-representative-sequences rep-seqs-deblur.qza --o-table table-deblur.qza --o-stats deblur-stats.qza) 2> DTime.txt


echo "************************** Genereating FeatureData summaries and Feature Tables *******************************************"
qiime feature-table summarize \
      --i-table table-deblur.qza \
      --o-visualization table.qzv \
qiime feature-table tabulate-seqs \
      --i-data rep-seqs-deblur.qza \
      --o-visualization rep-seqs.qzv

mkdir final

qiime tools export table-deblur.qza --output-dir final
biom convert -i final/feature-table.biom --to-tsv -o table.tsv

echo "*************************** Done running forward read deblur pipeline ***************************************************"
