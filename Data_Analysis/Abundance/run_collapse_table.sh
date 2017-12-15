#!/bin/bash

for i in *.tsv; do name=${i/.*/}; ./collapse_table.py -t $i -o $name"_col.tsv" -p True; done
