#!/usr/bin/python




import argparse
import re

def main():

    parser = argparse.ArgumentParser(
        description="Gets seq I.ds that they match",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-b", "--blast", metavar="Blast_file", type=str, help="Path to blastn result file", required=True)
    parser.add_argument("-f", "--fasta", metavar="Fasta_file", type=str, help="Path to fasta file to read in", required=True)
    args = parser.parse_args()

    #get Study, Pipeline and Filter Name
    heads=args.blast.split('/')
    Study=heads[int(len(heads)-2)]
    fields = re.split("([HML])",heads[int(len(heads)-1)])
    pipe = fields[0]
    filt = fields[1]+fields[2].split(".")[0]
    Fasta_Seqs=read_fasta(args.fasta)
    AsvtoMatch={}


    #have to get every 6th line and make sure its over 97% match

    with open(args.blast, "r") as blast:
        lc = 0

        for line in blast:
            if lc == 3 and "Fields:" not in line:
                lc = 0
            elif lc <= 4:
                lc += 1
            else:
                lc = 0
                #read in all the values from the line
                name,match,identity,align,mismatch,gap,qstart,qend,sstart,send,ev,bit = line.split("\t")
                print(name)
                if identity == "100.000":
                    if int(len(Fasta_Seqs[name])) == int(align):
                        AsvtoMatch[name]=match
                    else:
                        score=float(align)/float(len(Fasta_Seqs[name]))
                        if score >= 0.97:
                            AsvtoMatch[name]=match
                elif float(identity) >= 97.0:

                    if int(len(Fasta_Seqs[name])) == int(align):
                        AsvtoMatch[name]=match
                    else:
                        Addmismatch = int(len(Fasta_Seqs[name]))-int(align)
                        newmismatch = Addmismatch + int(mismatch)
                        score = 1.0-(float(newmismatch))/float(len(Fasta_Seqs[name]))
                        if float(score) >= 0.97:
                            AsvtoMatch[name]=match
    with open(Study+"_"+pipe+"_"+filt+".matchs", "w") as out:
        for name in AsvtoMatch:
            out.write(name+"\t"+AsvtoMatch[name]+"\n")
                        
                    
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
