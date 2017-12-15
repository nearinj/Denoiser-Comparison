#!/bin/bash

source activate qiime2-2017.9

echo "******************************** Reading in Config File *******************************************************"
echo "Reading in " $1

readarray -t vars <$1

if [ ${vars[0]} -eq 1 ]; then
    echo "Removing old data"
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/deblurF
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/filtered_fastqs/renamed
fi


cd ../../../DenoiseCompare_Out/${vars[1]}/$2


echo "******************************* Creating links to import into qiime2 artifact ********************************************"

cd filtered_fastqs

mkdir renamed

for i in *R1*fastq; do name=${i/_*/}; cp  $PWD/$i $PWD/renamed/$name"_X_L001_R1_001.fastq"; done
cd renamed

gzip *


cd ../../

echo "***************************** Importing filtered files into qiime2 artifact **********************************************"

mkdir deblurF
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path filtered_fastqs/renamed --source-format CasavaOneEightSingleLanePerSampleDirFmt --output-path deblurF/filtered.qza

cd deblurF
if [ ${vars[1]} == "Its-Mock9" ]; then
    echo "**************** Running deblur on ITS reads ***********************************************"
    (time qiime deblur denoise-other --p-jobs-to-start ${vars[4]} --i-demultiplexed-seqs filtered.qza --p-trim-length ${vars[2]} --o-representative-sequences rep-seqs-deblur.qza --o-table table-deblur.qza --o-stats deblur-stats.qza --output-dir extra --i-reference-seqs ../../../DB/ITs_DB.qza) 2>DTime.txt
else
    echo "******************** Running Deblur on forward reads ********************************************"
    (time qiime deblur denoise-16S --p-jobs-to-start ${vars[4]} --i-demultiplexed-seqs filtered.qza --p-trim-length ${vars[2]} --o-representative-sequences rep-seqs-deblur.qza --o-table table-deblur.qza --o-stats deblur-stats.qza --output-dir extra) 2>DTime.txt
fi

echo "************************** Genereating FeatureData summaries and Feature Tables *******************************************"
qiime feature-table summarize \
      --i-table table-deblur.qza \
      --o-visualization table.qzv \
qiime feature-table tabulate-seqs \
      --i-data rep-seqs-deblur.qza \
      --o-visualization rep-seqs.qzv
mkdir final

qiime tools export table-deblur.qza --out-dir final
biom convert -i final/feature-table.biom --to-tsv -o seqs.tsv

echo "*************************** Done running forward read deblur pipeline ***************************************************"
