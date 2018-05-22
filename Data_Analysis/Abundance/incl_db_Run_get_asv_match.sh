#!/bin/bash

#HMP DADA
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/high/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaHigh.out -b2 ../blast_nonmatch/HMP_Dada_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/low/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaLow.out -b2 ../blast_nonmatch/HMP_Dada_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/med/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaMed.out -b2 ../blast_nonmatch/HMP_Dada_Med_nonmatch.out

#HMP Deblur
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/high/deblurP/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurHigh.out -b2 ../blast_nonmatch/HMP_Deblur_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/low/deblurP/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurLow.out -b2 ../blast_nonmatch/HMP_Deblur_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/med/deblurP/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurMed.out -b2 ../blast_nonmatch/HMP_Deblur_Med_nonmatch.out

# HMP Unoise
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/high/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseHigh.out -b2 ../blast_nonmatch/HMP_Unoise_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/low/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseLow.out -b2 ../blast_nonmatch/HMP_Unoise_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/med/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseMed.out -b2 ../blast_nonmatch/HMP_Unoise_Med_nonmatch.out


# Mock-12 DADA
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/high/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaHigh.out -b2 ../blast_nonmatch/Mock-12_Dada_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/low/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaLow.out -b2 ../blast_nonmatch/Mock-12_Dada_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/med/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaMed.out -b2 ../blast_nonmatch/Mock-12_Dada_Med_nonmatch.out

# Mock-12 Deblur
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/high/deblurF/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurHigh.out -b2 ../blast_nonmatch/Mock-12_Deblur_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/low/deblurF/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurLow.out -b2 ../blast_nonmatch/Mock-12_Deblur_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/med/deblurF/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurMed.out  -b2 ../blast_nonmatch/Mock-12_Deblur_Med_nonmatch.out

# Mock-12 Unoise
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/high/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseHigh.out -b2 ../blast_nonmatch/Mock-12_Unoise_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/low/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseLow.out -b2 ../blast_nonmatch/Mock-12_Unoise_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/med/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseMed.out -b2 ../blast_nonmatch/Mock-12_Unoise_Med_nonmatch.out

# Mock-9 DADA
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/high/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaHigh.out -b2 ../blast_nonmatch/Mock-9_Dada_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/low/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaLow.out -b2 ../blast_nonmatch/Mock-9_Dada_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/med/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaMed.out -b2 ../blast_nonmatch/Mock-9_Dada_Med_nonmatch.out

#Mock-9 Deblur
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/high/deblurF/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurHigh.out -b2 ../blast_nonmatch/Mock-9_Deblur_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/low/deblurF/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurLow.out -b2 ../blast_nonmatch/Mock-9_Deblur_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/med/deblurF/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurMed.out -b2 ../blast_nonmatch/Mock-9_Deblur_Med_nonmatch.out

# Mock-9 Unoise
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/high/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseHigh.out -b2 ../blast_nonmatch/Mock-9_Unoise_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/low/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseLow.out -b2 ../blast_nonmatch/Mock-9_Unoise_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/med/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseMed.out -b2 ../blast_nonmatch/Mock-9_Unoise_Med_nonmatch.out

# Zymock Dada
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/high/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaHigh.out -b2 ../blast_nonmatch/Zymock_Dada_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/low/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaLow.out -b2 ../blast_nonmatch/Zymock_Dada_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/med/dada2/seqtab.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaMed.out -b2 ../blast_nonmatch/Zymock_Dada_Med_nonmatch.out

#Zymock Deblur
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/high/deblurP/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurHigh.out -b2 ../blast_nonmatch/Zymock_Deblur_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/low/deblurP/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurLow.out -b2 ../blast_nonmatch/Zymock_Deblur_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/med/deblurP/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurMed.out -b2 ../blast_nonmatch/Zymock_Deblur_Med_nonmatch.out

# Zymock Unoise
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/high/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseHigh.out -b2 ../blast_nonmatch/Zymock_Unoise_High_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/low/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseLow.out -b2 ../blast_nonmatch/Zymock_Unoise_Low_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/med/Unoise/zotus.fa -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseMed.out -b2 ../blast_nonmatch/Zymock_Unoise_Med_nonmatch.out

# OpenRef

./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Zymock/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Zymock/OpenMed.out -b2 ../blast_nonmatch/Zymock_Open_Med_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/HMPV4-V5/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/OpenMed.out -b2 ../blast_nonmatch/HMP_Open_Med_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Its-Mock9/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-9/OpenMed.out -b2 ../blast_nonmatch/Mock-9_Open_Med_nonmatch.out
./incl_db_get_asv_match.py -f ../../../../DenoiseCompare_Out/Mock-12/med/OpenRef_OTU/final/dna-sequences.fasta -b ../../../../DenoiseCompare_Out/Blast_Results/Mock-12/OpenMed.out -b2 ../blast_nonmatch/Mock-12_Open_Med_nonmatch.out


















