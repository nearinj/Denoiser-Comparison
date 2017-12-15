#!/bin/bash


#Use config file for Trim
source activate dada2

readarray -t vars <$1


echo "********************************************* Creating Filtered 5k ********************************************"

cd ../../../DenoiseCompare_Out/${vars[2]}/${vars[4]}


mkdir Filt5k

for i in filtered_fastqs/*; do
    name=${i#*/};
    head -n 20000 $i > Filt5k/$name;
done



echo "*************************************** Creating Filtered 10k *************************************************"

mkdir Filt10k

for i in filtered_fastqs/*; do
    name=${i#*/};
    head -n 40000 $i > Filt10k/$name;
done



mkdir Filt20k

for i in filtered_fastqs/*; do
    name=${i#*/};
    head -n 80000 $i > Filt20k/$name;
done


mkdir Filt30k

for i in filtered_fastqs/*; do
    name=${i#*/};
    head -n 120000 $i > Filt30k/$name;
done



mkdir Filt40k

for i in filtered_fastqs/*; do
    name=${i#*/};
    head -n 160000 $i > Filt40k/$name;
done


