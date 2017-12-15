#!/usr/bin/python3

#Renames fastq sample names
import argparse
import os


def main():

    parser = argparse.ArgumentParser(
        description="Rename Fastq's to be used in Usearch Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--fastq", metavar="FASTQ", type=str, help="Path to FastQ file", required=True)
    parser.add_argument("-s", "--sample", metavar="Sample_Name", type=str, help="Name of sample", required=True)
    parser.add_argument("-o", "--output", metavar="Output_Directory", type=str, help="Name of Output Directory", required=True)
    parser.add_argument("-r", "--rev", metavar="Reverse = true", type=str, help="set if doing reverse reads")
    args = parser.parse_args()

    #check if output exists
    if os.path.exists(args.output):
        print ("exists")
    else:
        os.makedirs(args.output)

    # initialize FastQ line counter
    fastq_lc = 4
    # initialize read number
    fastq_sn = 1
    with open(args.fastq, 'r') as file:
        data = file.readlines()

    i = 0
    while i != len(data):
        if i%4 == 0:
            data[i]='@' + args.sample + '.' + str(fastq_sn) + '\n'
            fastq_sn += 1
        i += 1
        
    if args.rev == "y":
        outfile=args.output + '/' + args.sample + "_R2_rename.fastq"
    else:
        outfile=args.output + '/' + args.sample + "_R1_rename.fastq"
            
    with open(outfile, 'w') as file:
        file.writelines( data )
        
    
if __name__ == '__main__':
        main()
            
    
    
        

    
    

    
