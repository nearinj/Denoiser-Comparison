#!/bin/bash

for i in *.tsv; do
    name=${i/.*/};
    echo $name
    ./collapse_table.py -t $i -o $name"_col.tsv" -p True;
done
