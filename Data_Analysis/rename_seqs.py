#!/usr/bin/python

import argparse
import re


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", metavar="fasta file", type=str, required=True)
    parser.add_argument("-o", "--out", metavar="out fastafile", type=str, required=True)
    args=parser.parse_args()



    with open(args.input, "r") as infile:
        with open (args.out, "w") as outfile:
            lc=1
            i=1;
            for line in infile:
                if lc == 1:
                    info= re.sub(">","",line)
                    infos=info.rstrip("\n")
                    outfile.write(infos+"\t")
                    lc += 1
                    i += 1
                elif lc == 2:
                    outfile.write(line)
                    lc =1;
                    
                    
if __name__ == '__main__':
    main()
    
