#!/bin/bash


#all the HMP Data
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/high/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/high/deblurP/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/high/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/med/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/med/deblurP/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/med/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/low/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/low/deblurP/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/HMPV4-V5/low/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseLow.out -o BlastResults.tsv


#all the Zymock Data
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/high/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Zymock/DadaHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/high/deblurP/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/high/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/med/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Zymock/DadaMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/med/deblurP/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/med/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/low/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Zymock/DadaLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/low/deblurP/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Zymock/low/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseLow.out -o BlastResults.tsv



#all the Mock-12 Data
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/high/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/high/deblurF/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/high/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/med/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/med/deblurF/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/med/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/low/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/low/deblurF/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Mock-12/low/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseLow.out -o BlastResults.tsv


#all the Mock-9 Data

./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/high/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/high/deblurF/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/high/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseHigh.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/med/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/med/deblurF/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/med/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseMed.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/low/dada2/seqtab.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/low/deblurF/final/dna-sequences.fasta -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurLow.out -o BlastResults.tsv
./get_blast_results.py -f ../../DenoiseCompare_Out/Its-Mock9/low/Unoise/zotus.fa -b ../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseLow.out -o BlastResults.tsv



