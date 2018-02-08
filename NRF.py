import pysam
import sys

if sys.argv.__len__() != 2:
    print("python2.7 NRF.py <input_file>")
else:
    fname_in = sys.argv[1]
    infile = pysam.AlignmentFile(fname_in, "rb")
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

    print(uni, tot, uni/tot)
