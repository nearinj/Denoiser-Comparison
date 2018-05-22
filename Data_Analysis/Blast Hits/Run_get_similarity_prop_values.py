#!/bin/bash


#all the HMP Data
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaHigh.out -o Sims/HMP_DadaHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/high/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurHigh.out -o Sims/HMP_DeblurHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseHigh.out -o Sims/HMP_UnoiseHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaMed.out -o Sims/HMP_DadaMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurMed.out -o Sims/HMP_DeblurMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseMed.out -o Sims/HMP_UnoiseMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaLow.out -o Sims/HMP_DadaLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/low/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurLow.out -o Sims/HMP_DeblurLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseLow.out -o Sims/HMP_UnoiseLow.tsv


#all the Zymock Data
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaHigh.out -o Sims/Zymock_DadaHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/high/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurHigh.out -o Sims/Zymock_DeblurHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseHigh.out -o Sims/Zymock_UnoiseHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaMed.out -o Sims/Zymock_DadaMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/med/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurMed.out -o Sims/Zymock_DeblurMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseMed.out -o Sims/Zymock_UnoiseMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaLow.out -o Sims/Zymock_DadaLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/low/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurLow.out -o Sims/Zymock_DeblurLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseLow.out -o Sims/Zymock_UnoiseLow.tsv



#all the Mock-12 Data
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaHigh.out -o Sims/Mock12_DadaHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/high/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurHigh.out -o Sims/Mock12_DeblurHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseHigh.out -o Sims/Mock12_UnoiseHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaMed.out -o Sims/Mock12_DadaMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/med/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurMed.out -o Sims/Mock12_DeblurMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseMed.out -o Sims/Mock12_UnoiseMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaLow.out -o Sims/Mock12_DadaLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/low/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurLow.out -o Sims/Mock12_DeblurLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseLow.out -o Sims/Mock12_UnoiseLow.tsv


#all the Mock-9 Data

./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaHigh.out -o Sims/Mock9_DadaHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/high/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurHigh.out -o Sims/Mock9_DeblurHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseHigh.out -o Sims/Mock9_UnoiseHigh.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaMed.out -o Sims/Mock9_DadaMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurMed.out -o Sims/Mock9_DeblurMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseMed.out -o Sims/Mock9_UnoiseMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaLow.out -o Sims/Mock9_DadaLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/low/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurLow.out -o Sims/Mock9_DeblurLow.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseLow.out -o Sims/Mock9_UnoiseLow.tsv



#get all the Open-ref results... need to do this.
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP//HMP/OpenMed.out -o Sims/HMP_OpenMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Zymock/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/OpenMed.out -o Sims/Zymock_OpenMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/OpenMed.out -o Sims/Mock9_OpenMed.tsv
./get_similarity_prop_values.py -f ../../../DenoiseCompare_Out/Mock-12/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/OpenMed.out -o Sims/Mock12_OpenMed.tsv

#done!
