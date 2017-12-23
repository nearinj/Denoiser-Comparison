#!/usr/bin/python

#reads in time directory that has the time file in it and extracts the time data while converting it to minutes

import argparse
import re
import os.path
from collections import deque

def main():
    parser=argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", metavar="Directory_of_filt", type=str, help="Path to filt Directory", required=True)
    parser.add_argument("-o", "--output", metavar="File to output", type=str, help="Path to filt Directory", required=True)




    
    args = parser.parse_args()



    #get Rarifcation number
    rare= re.split("([0-9]+k)",args.directory)
    number= rare[1]
    

    if os.path.isfile(args.output) is False:
        with open(args.output, "w") as infile:
            infile.write("Filt\tPipe\tTimeCMD\tReal\tUser\tSys\n")

    #Get all the dictionaries together

    
    DadaCT = {}
    DadaCT = getTimeInfo(DadaCT, "Dada", "timeCT", args.directory)

    DadaI = {}
    DadaI = getTimeInfo(DadaI, "Dada", "timeI", args.directory)

    deblurD = {}
    deblurD = getTimeInfo(deblurD, "Deblur", "DTime", args.directory)

    UnoiseD = {}
    UnoiseD = getTimeInfo(UnoiseD, "Unoise", "Dtime", args.directory)

    UnoiseUni = {}
    UnoiseUni = getTimeInfo(UnoiseUni, "Unoise", "Unitime", args.directory)

    All_Dic = [DadaCT, DadaI, deblurD, UnoiseD, UnoiseUni]
    
    with open(args.output, "a") as infile:
        for Dic in All_Dic:
            infile.write(number+"\t"+Dic["pipe"]+"\t"+Dic["TimeCMD"]+"\t"+str(Dic["real"])+"\t"+str(Dic["user"])+"\t"+str(Dic["sys"])+"\n")
            
            



def getTimeInfo(Dic, Pipe, TimeCMD, path):
    Dic["pipe"] = Pipe
    Dic["TimeCMD"] = TimeCMD

    if Pipe == "Dada":
        Address1 = "Dadaout"
    elif Pipe == "Deblur":
        Address1 = "deblurP"
    elif Pipe == "Unoise":
        Address1 = "Unoise_out"

    with open (path+"/"+Address1+"/"+TimeCMD+".txt", "r") as infile:
        last3 = deque(infile, 3)
        realtime = last3[0].split("\t")[1]
        Dic["real"] = convertTime(realtime)

        usertime = last3[1].split("\t")[1]
        Dic["user"] = convertTime(usertime)

        Systime = last3[2].split("\t")[1]
        Dic["sys"] = convertTime(Systime)

    return Dic


def convertTime(time):
    minutes=time.split("m")
    StoM=float(minutes[0])*60
    seconds=minutes[1].replace("s", "")
    FixedTime=float(seconds)+StoM

    return FixedTime


if __name__ == '__main__':
    main()
