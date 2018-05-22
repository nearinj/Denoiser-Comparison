#!/bin/bash


for i in *.fasta; do
    name=${i/.*/};
    echo $name;
    if [[ $name == "Mock-9"* ]]; then
       blastn -db ../../blast_db/Unite/sh_general_release_dynamic_s_10.10.2017_dev.fasta -query $i -out blast_nonmatch/$name.out -max_target_seqs 1 -outfmt 7
       echo "blasting against Silva"
    else
	blastn -db ../../blast_db/silva/SILVA_128_SSURef_Nr99_tax_silva.fasta -query $i -out blast_nonmatch/$name.out -max_target_seqs 1 -outfmt 7;
    fi
done
