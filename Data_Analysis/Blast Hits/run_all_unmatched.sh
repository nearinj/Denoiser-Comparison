#!/bin/bash



for i in *.fasta; do name=${i/.*/}; echo $name; ./get_unmatched_result.py -f $i -b blast_nonmatch/$name.out -o Silva_Blast_Results.tsv; done 
