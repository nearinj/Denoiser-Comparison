#!/bin/bash



source activate qiime2-2017.12

echo "**************************************************Running Denovo-OTU picking using Vsearch QIIME2-plugin*************************************"


echo "Reading in " $1

readarray -t vars <$1


if [ "${vars[0]}" == "1" ]; then
    echo "Removing old data"
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/OpenRef_OTU/
fi


cd ../../../DenoiseCompare_Out/${vars[1]}/$2

echo "../../DenoiseCompare_Out/"${vars[1]}"/"$2




echo "*************************************** Run this script after the Deblur script or paired files will be missing!!! *****************************************"

mkdir OpenRef_OTU
echo $PWD

echo "******************** converting to Qiime2 artifact ******************************"

qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path filtered_fastqs/renamed --output-path OpenRef_OTU/filtered.qza --source-format CasavaOneEightSingleLanePerSampleDirFmt 

cd OpenRef_OTU

echo "**************************** Dereplicating sequences ****************************************"

(/usr/bin/time -v qiime vsearch dereplicate-sequences --i-sequences filtered.qza --o-dereplicated-table derep-table.qza --o-dereplicated-sequences derep-seqs.qza) 2> DerepTime.txt

if [ ${vars[1]} == "Its-Mock9" ]; then

    echo "********************************* removing ITS Chimeras ***********************************************"

    (/usr/bin/time -v qiime vsearch uchime-ref --i-sequences derep-seqs.qza --i-table derep-table.qza --i-reference-sequences ~/projects/DenoiseCompare/blast_db/Unite/chimera_db/97_Unite_db.qza --o-chimeras chimeras.qza --o-nonchimeras non-chim-derep-seqs.qza --o-stats chimera-stats.qza --p-threads ${vars[3]}) 2> ChimeraTime.txt

    qiime feature-table filter-features --m-metadata-file non-chim-derep-seqs.qza --i-table derep-table.qza --o-filtered-table non-chim-table.qza

    echo "********************************** Rumning ITS open-ref picking ***************************************"
    (/usr/bin/time -v qiime vsearch cluster-features-open-reference --i-table non-chim-table.qza --i-sequences non-chim-derep-seqs.qza --o-clustered-table open-ref-table-97.qza --o-clustered-sequences rep-seqs-open-ref-97.qza --p-perc-identity 0.97 --p-threads ${vars[3]} --i-reference-sequences ~/projects/DenoiseCompare/blast_db/Unite/chimera_db/97_Unite_db.qza --o-new-reference-sequences new_db.qza) 2> ClustTime.txt

elif [ ${vars[1]} != "Its-Mock9" ]; then
    echo "********************************************************************** Remove Chimeras 16S****************************************************"
    (/usr/bin/time -v qiime vsearch uchime-ref --i-sequences derep-seqs.qza --i-table derep-table.qza --i-reference-sequences ~/projects/DenoiseCompare/blast_db/GreenGenes13_8_97/gg_13_8_otus/rep_set/97_otus.qza --o-chimeras chimeras.qza --o-nonchimeras non-chim-derep-seqs.qza --o-stats chimera-stats.qza --p-threads ${vars[3]}) 2> ChimeraTime.txt

    qiime feature-table filter-features --m-metadata-file non-chim-derep-seqs.qza --i-table derep-table.qza --o-filtered-table non-chim-table.qza
    echo "************************************ Running Open Ref OTU picking ****************************************"
    (/usr/bin/time -v qiime vsearch cluster-features-open-reference  --i-table non-chim-table.qza --i-sequences non-chim-derep-seqs.qza --o-clustered-table open-ref-table-97.qza --o-clustered-sequences rep-seqs-open-ref-97.qza --p-perc-identity 0.97 --p-threads ${vars[3]} --i-reference-sequences ~/projects/DenoiseCompare/blast_db/Unite/chimera_db/97_Unite_db.qza --o-new-reference-sequences new_db.qza ) 2> ClustTime.txt


fi

echo "********************************************** removing singletons ***************************************"
qiime feature-table filter-features --i-table open-ref-table-97.qza --p-min-frequency 2 --output-dir remove_one
qiime feature-table filter-seqs --i-data rep-seqs-open-ref-97.qza --i-table remove_one/filtered_table.qza --o-filtered-data remove_one/rep-seqs_filtered.qza

echo "******************************************** Generating FeatureData summaries and Feature tables ***********************************"
qiime feature-table summarize --i-table open-ref-table-97.qza --o-visualization open-ref-table-97.qzv

mkdir final

qiime tools export open-ref-table-97.qza --output-dir final


echo "*************************************** Done running Open Ref OTU picking ***********************************************************"

