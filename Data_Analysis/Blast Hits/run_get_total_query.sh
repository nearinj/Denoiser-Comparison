#!/bin/bash


#all the HMP Data
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DadaLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/DeblurLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/UnoiseLow.out -o Total_query.tsv


#all the Zymock Data
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DadaLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/DeblurLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/UnoiseLow.out -o Total_query.tsv



#all the Mock-12 Data
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DadaLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/DeblurLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/UnoiseLow.out -o Total_query.tsv


#all the Mock-9 Data

./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseHigh.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DadaLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/DeblurLow.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/UnoiseLow.out -o Total_query.tsv



#get all the Open-ref results... need to do this.
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/HMP/HMP/OpenMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Zymock/OpenMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-9/OpenMed.out -o Total_query.tsv
./get_total_query_seq.py -b ../../../DenoiseCompare_Out/Blast_Results/Mock-12/OpenMed.out -o Total_query.tsv

#done!
