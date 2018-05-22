#!/usr/bin/python

import argparse
import re
import os.path

def main():

    parser = argparse.ArgumentParser(
        description="Convert results of blastn into an R compatible tsv file",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--blast", metavar="Blast_file", type=str, help="Path to blastn result file", required=True)
    parser.add_argument("-f", "--fasta", metavar="fasta_file", type=str, help="Path to fasta file containning query sequences", required=True)
    parser.add_argument("-o", "--output", metavar="Output_file", type=str, help="Name of output file", required=True)
    args = parser.parse_args()


    #get Pipeline name and Filter Name
    heads=args.blast.split('/')
    fields = re.split("([HML])",heads[int(len(heads)-1)])
    pipe = fields[0]
    filt = fields[1]+fields[2].split(".")[0]
    sample = heads[int(len(heads)-2)]
    Pnummatch = 0
    APnummatch = 0
    NonmatchSeq = {}
    prevLine = ""
    prev2Line = ""
    test_miss= 0
    
    Fasta_Seqs=read_fasta(args.fasta)
    with open(args.output, "a") as out:
        with open(args.blast, "r") as blast:
            lc=0

            for line in blast:
                if lc == 3 and "Fields:" not in line:
                    lc = 0
                    if "0 hits" in line:
                        missed_seq_name=prev2Line.split(" ")[2].rstrip()
                        NonmatchSeq[missed_seq_name] = Fasta_Seqs[missed_seq_name]
                    
                   # values = prev2Line.split(" ")
                   # name = values[int(len(values)-1)]
                   # NonmatchSeq[name] = Fasta_Seqs[name]
                elif lc <= 4:
                    lc += 1
                else:
                    hits_found=prevLine.split(" ")[1]
                    lc = -int(hits_found)+1
                    
                    #read in all the values from the line
                    name,match,identity,align,mismatch,gap,qstart,qend,sstart,send,ev,bit = line.split("\t")
                    #check if identity is 100% and if so check if the alignment length = the length of the fasta sequence
                    if float(identity) == 100.00:
                        if int(len(Fasta_Seqs[name])) == int(align):
                            Pnummatch += 1
                        else:
                            score=float(float(qend)-float(qstart)+1)/float(len(Fasta_Seqs[name]))
                            if score >= 0.97:
                                
                                APnummatch += 1
                            else:
                                NonmatchSeq[name] = Fasta_Seqs[name]
                                test_miss += 1
                    elif float(identity) >= 97.0:
                        
                        if int(len(Fasta_Seqs[name])) == int(align):
                            
                            APnummatch += 1
                        else:
                            
                            not_aligned = float(len(Fasta_Seqs[name])) - (float(qend)-float(qstart)+1.0)
                            score = 1.0-(float(mismatch)+float(not_aligned)+float(gap))/float(len(Fasta_Seqs[name]))
                            if float(score) >= 0.97:
                                
                                APnummatch += 1
                            else:
                                NonmatchSeq[name] = Fasta_Seqs[name]
                                test_miss +=1
                    else:
                        NonmatchSeq[name] = Fasta_Seqs[name]
                        test_miss +=1
                    print(name)
                prev2Line=prevLine
                prevLine=line
                                
        out.write("\n"+sample+"\t"+pipe+"\t"+filt+"\t"+str(Pnummatch)+"\t"+str(APnummatch))
                            
    with open(sample+"_"+pipe+"_"+filt+"_"+"nonmatch.fasta", "w") as nonmatch:
        for name in NonmatchSeq:
            nonmatch.write(">"+name+"\n "+NonmatchSeq[name]+"\n")
                
    print(len(NonmatchSeq))
    print(test_miss)

    #function written by Gavin Douglas
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




