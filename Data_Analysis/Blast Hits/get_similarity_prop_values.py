#!/usr/bin/python


import argparse
import re
import os.path

def main():

    parser=argparse.ArgumentParser(description="Takes a blast file and a fasta file and outputs a list of identity precentages taking into account the length of the sequences in the fasta")
    parser.add_argument("-b", "--blast", metavar="Blast_file", type=str, required=True)
    parser.add_argument("-f", "--fasta", metavar="fasta_file", type=str, required=True)
    parser.add_argument("-o", "--output", metavar="output_file", type=str, required=True)

    args = parser.parse_args()




    #read fasta file
    Fasta_seqs=read_fasta(args.fasta)
    #make dictionary to store name and identity % match
    identity_dic={}
    total_query_seqs=0
    with open(args.output, "w") as out:
        with open(args.blast, "r") as blast:
            lc = 0

            for line in blast:
                if "# BLAST processed " in line:
                    total_query_seqs=line.split(" ")[3]
                if lc == 3 and "Fields:" not in line:
                    lc = 0
                elif lc <= 4:
                    lc += 1
                else:
                    lc = 0
                    #read in all the values from the line
                    name,match,identity,align,mismatch,gap,qstart,qend,sstart,send,ev,bit = line.split("\t")
                    #check if identity is 100% and if so check if alignment length = the length of the fasta sequence
                    if identity == "100.000":
                        if(int(len(Fasta_seqs[name]))) == int(align):
                           identity_dic[name]=identity
                        else:
                           score=float(float(qend)-float(qstart)+1.0)/float(len(Fasta_seqs[name]))
                           identity_dic[name]=score*100
                    else:
                        if(int(len(Fasta_seqs[name]))) == int(align):
                           identity_dic[name]=identity
                        else:
                            not_aligned= float(len(Fasta_seqs[name]))-float(float(qend)-float(qstart)+1.0)
                            score = 1.0-(float(not_aligned)+float(gap)+float(mismatch))/float(len(Fasta_seqs[name]))
                            identity_dic[name]=score*100
                            
        for i in identity_dic:
            out.write(i+"\t"+str(identity_dic[i])+"\n")
        for i in range(int(total_query_seqs) - int(len(identity_dic))):
            out.write("no_hit"+str(i)+"\t"+"0\n")
    if(int(total_query_seqs) != int(len(identity_dic))):
        print("some sequences had 0 hits ! in"+args.blast)
        
def read_fasta(filename):


    seq = {}
    name = None

    with open(filename, "r") as fasta:

        for line in fasta:

            if line[0] == ">":
                name = line.split()[0][1:]
                seq[name] = ""

            else:

                line = line.rstrip("\r\n")
                seq[name] += line
    return seq


if __name__ == '__main__':
    main()
                                                                                
