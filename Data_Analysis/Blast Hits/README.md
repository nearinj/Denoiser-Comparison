### blast_nonmatch.sh
Takes all of the fasta files that contain sequences that did not match with the expected sequences at 97% or greater and runs them on a blast search against the SILVA 16S rRNA gene database. 
### get_blast_results.py
Takes two parameters a blast file in tabular format output and the corresponding fasta file that was used to generate it. It then outputs a tab delimited table with the amount of 100% and 97% expected hits from that blast search as well fasta files that contain all of the sequences that did not match with the expected sequeces at 97% or greater for each pipeline and filter stringency.
### get_unmatched_result.py
Takes the resulting blast file from blast_nonmatch.sh and the fasta file that was used to generate it as parameters. It then outputs a tab delimited table with the amount of 100% and 97% or greater matchs.
### run_all_unmatched.sh
Runs get_unmatched_result.py on all of the 
### run_get_blast_results.sh
Runs the get_blast_results.py on all of the data generated from the mock communities
