#!/usr/bin/python


import argparse
import re
import os.path

def main():

    parser= argparse.ArgumentParser(description="Get total amount of sequences blasted")
    parser.add_argument("-b", "--blast", metavar="Blast_file", type=str,  help="Path to blast file you would like to know the total amount of querys for", required=True)
    parser.add_argument("-o", "--output", metavar="Blast_file", type=str, help="output for the total amount")

    args=parser.parse_args()
    
    #get pipeline name and Filter name
    heads=args.blast.split("/")
    fields=re.split("([HML])",heads[int(len(heads)-1)])
    pipe=fields[0]
    filt=fields[1]+fields[2].split(".")[0]
    sample= heads[int(len(heads)-2)]






    fileHandle = open(args.blast, "r")
    linelist = fileHandle.readlines()
    fileHandle.close()
    
    total_query = linelist[len(linelist)-1].split(" ")[3]
    
    with open(args.output, "a") as out:
        out.write(sample+"\t"+pipe+"\t"+filt+"\t"+total_query+"\n")



if __name__ == "__main__":
    main()
    
