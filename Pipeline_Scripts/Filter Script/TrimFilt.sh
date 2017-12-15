#!/bin/bash

echo "Reading in "$1
readarray -t vars <$1

source activate dada2

echo "**************************************************************** Filtering and Trimming *********************************************"

if [ ${vars[1]} -eq 1 ]; then
    echo "removing old information"
    rm -r ../../../DenoiseCompare_Out/${vars[2]}/${vars[4]}
fi

echo "making HMP directory"
mkdir ../../../DenoiseCompare_Out/${vars[2]}

cd ../../../DenoiseCompare_Out/${vars[2]}

echo "making " ${vars[4]} " directory"

mkdir ${vars[4]}

cd ${vars[4]}

mkdir primer_trimmed_fastqs

echo ${vars[5]}
echo ${vars[6]}

if [[ ${vars[0]} != "skipcut" ]]; then
#export the array to be used in GNU parallel
export vars
parallel --env vars --link --jobs ${vars[3]} --eta \
	   'cutadapt \
    --pair-filter any \
    --no-indels \
    --discard-untrimmed \
    -g '${vars[5]}' \
    -G '${vars[6]}' \
    -o primer_trimmed_fastqs/{1/}.gz \
    -p primer_trimmed_fastqs/{2/}.gz \
    {1} {2} \
    > primer_trimmed_fastqs/{1/}_cutadapt_log.txt' \
	     ::: ~/Mock_Communities/MockData/${vars[2]}/*_R1*.fastq ::: ~/Mock_Communities/MockData/${vars[2]}/*_R2*.fastq

echo "Making cut logs"

#put cut logs together
parse_cutadapt_logs.py -i primer_trimmed_fastqs/*log.txt
mkdir logs
mv primer_trimmed_fastqs/*log.txt logs

fi
echo "Picking filtering settings"

#set the filter based on file input
if [[ ${vars[4]} == "low" ]]; then
    ee=5
fi
if [[ ${vars[4]} == "med" ]]; then
    ee=3
fi
if [[ ${vars[4]} == "high" ]]; then
    ee=1
fi
echo "Picked "${vars[4]}

if [[ ${vars[0]} == "skipcut" ]]; then
    
    echo "Skipping Primer Trimming"
    dada2_filter.R -f ~/Mock_Communities/MockData/${vars[2]} --truncLen ${vars[7]} --maxN 0 --maxEE $ee --threads ${vars[3]} --f_match _R1.*fastq -s
    echo "dada2_filter.R -f ~/Mock_Communities/MockData/"${vars[2]}" --truncLen "${vars[7]}","${vars[8]}" --maxN 0 --maxEE "$ee","$ee" --threads "${vars[3]}
else

    echo "*************************************** Filtering Reads **********************************************"

    dada2_filter.R -f primer_trimmed_fastqs --truncLen ${vars[7]},${vars[8]} --maxN 0 --maxEE $ee,$ee --truncQ 2 --threads ${vars[3]} --f_match _R1.fastq.gz --r_match _R2.fastq.gz
fi

cd filtered_fastqs
gunzip *

echo "***************************************** Finished Trimming and Filtering reads ******************************"
