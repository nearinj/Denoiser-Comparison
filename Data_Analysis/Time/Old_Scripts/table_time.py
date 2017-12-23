#!/usr/bin/python


#reads in a time file and puts it into a table.
import sys

lc = 0

info = []

print("dataset	filter	pipeline	time_cmd	real	user	sys")

with open(sys.argv[1], "r") as infile:
	
	for line in infile:
	
		if "==>" == line[0:3]:
			
			if lc == 0:
				lc += 1
			else:
		                print("\t".join(info))
				info = []

			line_info = line.split(" ")[1]
			line_info = line_info.replace("./", "")
			line_info = line_info.replace(".txt", "")
			info = line_info.split("/")
			
		elif len(line) > 1:
			line = line.rstrip("\n")
			line_split = line.split("\t")
                        minutes=line_split[1].split("m")
                        StoM=float(minutes[0])*60
                        seconds=minutes[1].replace("s","")
                        time=float(seconds)+StoM
                        
			info += [str(time)]

print("\t".join(info))
