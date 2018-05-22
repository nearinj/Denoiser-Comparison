#!/bin/bash

echo "***************************** Running deblur without postive filt **************************************"
echo "Reading in "$1

readarray -t vars <$1

if [ ${vars[1]} == "1" ]; then
    echo "removing old data"
    rm -rf ../../../DenoiseCompare_Out/${vars[1]}/$2/Deblur_NoPos
fi


cd ../../../DenoiseCompare_Out/${vars[1]}/$2
echo $PWD
mkdir Deblur_NoPos
echo "************************* Run this after running the main deblur script **********************************"

cd Deblur_NoPos

(/usr/bin/time -v deblur workflow --seqs-fp ../stitched_reads/renamedpaired --output-dir output -t 350) 2> DeblurTime.txt

echo "*************************** Done running deblur without postive filtering ******************************"
