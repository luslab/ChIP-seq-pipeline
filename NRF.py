import pysam
import sys

infile = pysam.AlignmentFile(snakemake.input[0], "rb")
unique = set()
tot = 0
uni = 0
for line in infile:
    string = str(line.reference_name) + "_" + str(line.reference_start)
    # print(string)
    if string in unique:
        tot += 1
    else:
        unique.add(string)
        uni += 1
        tot += 1

with open(snakemake.output[0], "w") as outfile:
	outfile.write(str({"Unique" : uni, "Total" : tot, "Ratio" : uni/tot}))
