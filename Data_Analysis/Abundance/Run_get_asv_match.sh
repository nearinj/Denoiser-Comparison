#!/bin/bash

#HMP DADA
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaMed.out

#HMP Deblur
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/high/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/low/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurMed.out

# HMP Unoise
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/HMPV4-V5/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseMed.out


# Mock-12 DADA
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaMed.out

# Mock-12 Deblur
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/high/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/low/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/med/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurMed.out

# Mock-12 Unoise
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Mock-12/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseMed.out

# Mock-9 DADA
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaMed.out

#Mock-9 Deblur
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/high/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/low/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/deblurF/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurMed.out

# Mock-9 Unoise
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Its-Mock9/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseMed.out

# Zymock Dada
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/high/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/low/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/med/dada2/seqtab.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaMed.out

#Zymock Deblur
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/high/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/low/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/med/deblurP/final/dna-sequences.fasta -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurMed.out

# Zymock Unoise
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/high/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseHigh.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/low/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseLow.out
./get_asv_match.py -f ../../../DenoiseCompare_Out/Zymock/med/Unoise/zotus.fa -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseMed.out



















