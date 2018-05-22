#!/bin/bash




for i in matchs/*; do studyname=${i/_*/}; study=${studyname##*/}; filtname=${i/.*/}; filt=${filtname##*_}; pipename=${i/_[HLM]*/}; pipe=${pipename##*_};
		      echo $i;
		      ./add_organism.py -b TSV_Bioms/$study/$pipe$filt".tsv" -t $i -o $study"_"$pipe"_"$filt"_Org.tsv";
done

