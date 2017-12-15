#!/usr/bin/python




import argparse
import re
from collections import deque

def main():

    parser=argparse.ArgumentParser()
    parser.add_argument("-t", "--table", metavar="TSV_Table", type=str, help="Path to table with organisms", required=True)
    parser.add_argument("-o", "--output", metavar="Output_Table", type=str, help="Path to output collapsed table", required=True)
    parser.add_argument("-p", "--percent", metavar="Output_Table in Percent", type=bool, help="Output the table as a percent as well",  required=False)
    args=parser.parse_args()
    
    line1=""
    line2=""
    dic = {}
    samplenum=""
    with open(args.table, "r") as in_file:
        lc = 0
        
        for line in in_file:
            if lc == 0:
                lc += 1
                line1=line
            elif lc == 1:
                lc += 1
                info=line.split("\t")
                samplenum=int(len(info))-2
                line2=line
            else:
                info=line.split("\t")
                otu=info[0]
                org=info[int(len(info))-1]
                temporg=re.split("([0-9]*)$",org)
                organism=temporg[0].rstrip("\n")
                pattern = re.compile("[A-Za-z]") 
                if not pattern.match(organism):
                    organism = "Non_Reference"
                #if organism is already in dictionary append the otu to its value
                if organism in dic:
                    otus=dic[organism]
                    otus.append(otu)
                    dic[organism]=otus
                    
                #if organism is not in dictionary make a list with the otu as first element and make it the org value
                else:
                    otus= [otu]
                    dic[organism]=otus
        #really bad solution n^2 time zzzz >.>
        with open(args.output, "w") as out_file:
            with open(args.table, "r") as in_file:
                #Print header file
                j = 1
                Samplenames = []
                while j <= samplenum:
                    Samplenames.append("Sample"+str(j))
                    j += 1
                out_file.write("Organism\tOtus")
                for sample in Samplenames:
                    out_file.write("\t"+sample)
                out_file.write("\n")
                #print info for each organism

                for organism in dic:
                    with open(args.table, "r") as in_file:
                        sampleValues= []
                        for line in in_file:
                            info=line.split("\t")
                            otu=info[0]
                            #if the otu is in the current organism add its sample values to the current sample values
                            if otu in dic[organism]:
                                i = 0
                                #get all the sample values
                                if len(sampleValues) == 0:
                                    while i < samplenum:
                                        sampleValues.append(info[i+1])
                                        i += 1
                                else:
                                    while i < len(sampleValues):
                                        tempval = 0
                                        tempval = float(sampleValues[i])
                                        tempval += float(info[i+1])
                                        sampleValues[i] = tempval
                                        i += 1
                    #write out the organism name
                    out_file.write(organism+"\t")
                    #write out all the otu's that match that organism
                    if organism =="No_Match":
                        out_file.write("X")
                    else:
                        for otu in dic[organism]:
                            out_file.write(otu+"-")
                    #write out the sample values for that organism
                    for SampleValue in sampleValues:
                        out_file.write("\t"+str(SampleValue))
                    out_file.write("\n")
    sampleTotals = []
    with open(args.table, "r") as in_file:
        lc = 0
        for line in in_file:
            if lc <= 1:
                lc += 1
            else:
                info=line.split("\t")
                i = 0
                if len(sampleTotals) == 0:
                    while i < samplenum:
                        sampleTotals.append(info[i+1])
                        i += 1
                else:
                    while i < samplenum:
                        tempval = 0
                        tempval = float(sampleTotals[i])
                        tempval += float(info[i+1])
                        sampleTotals[i] = tempval
                        i += 1
    with open(args.output, "a") as out_file:
        out_file.write("Total\tNA")
        for val in sampleTotals:
            out_file.write("\t"+str(val))
        out_file.write("\n")
    #Output table in percents with the expected percent as well
    if (args.percent == True):
        print("making percent table")
        with open(args.output, "r") as in_file:
            print("opened in_file")
            #get the totals for each sample
            sampletotals = []
            lastline = deque(in_file, 1)[0]
            infos = lastline.split("\t")
            i = 0
            filename = args.table.split("_")
            studyname = filename[0]
            while i < samplenum:
                sampletotals.append(infos[i+2].rstrip("\n"))
                i += 1
            print(sampletotals)
        with open(args.output, "r") as in_file:
            with open("per_"+args.output, "w") as out_file:
                lc = 0
                for line in in_file:
                    if lc == 0:    
                        out_file.write(line.rstrip("\n")+"\n")
                        lc += 1
                    else:
                        print(line)
                        info=line.split("\t")
                        #get the percent val for the samples
                        sampleval = []
                        j = 1
                        print(samplenum)
                        while j <= samplenum:
                            print("info[j+1]: "+info[j+1])
                            sampleval.append(info[j+1])
                            j += 1
                        print(sampleval)
                        percentval = []
                        k = 0
                        while k < int(len(sampleval)):
                            percent = float(sampleval[k])/float(sampletotals[k])
                            percentval.append(percent)
                            k += 1
                        #find the expected based on study name
                        out_file.write(info[0]+"\t"+info[1])
                        for val in percentval:
                            out_file.write("\t"+str(val))
                        out_file.write("\n")
                    
if __name__ == '__main__':
    main()
    

