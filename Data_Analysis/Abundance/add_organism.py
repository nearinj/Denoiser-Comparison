#!/usr/bin/python


import argparse



def main():


    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--biom", metavar="biom_in_tsv", type=str, help="Path to TSV format biom table", required=True)
    parser.add_argument("-t", "--table", metavar="Table_matchs", type=str, help="Table with otu and species matches", required=True)
    parser.add_argument("-o", "--output", metavar="Output_Table", type=str, help="Path to output table", required=True)
    args = parser.parse_args()

    #get study, pipe, and filt name
    fields=args.table.split("_")
    study=fields[0]
    pipe=fields[1]
    filt=fields[2]

    lines= {}
    
    #Read in table file
    with open(args.table, "r") as table:
        for line in table:
            info=line.split("\t")
            otu=info[0]
            species=info[1].rstrip("\n")

            lines[otu]=species

    

    with open(args.output, "w") as out_file:
        with open(args.biom, "r") as in_file:
            lc=0
            for line in in_file:
                if lc == 0:
                    out_file.write(line)
                    lc += 1
                elif lc == 1:
                    newline=line.rstrip("\n")
                    out_file.write(newline+"\tOrganism\n")
                    lc += 1
                else:
                    info=line.split("\t")
                    if info[0] in lines:
                        out_file.write(line.rstrip('\n') + "\t" + lines[info[0]] + '\n')
                    else:
                        out_file.write(line)



if __name__ == '__main__':
    main()
    
