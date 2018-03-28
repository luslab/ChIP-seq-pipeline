import pysam
import sys

if sys.argv.__len__() != 2:
    print("python2.7 cc.py <input_file>")
else:
    fname_in = sys.argv[1]
    infile = pysam.AlignmentFile(fname_in, "rb")
    pos = 0
    neg = 0
    for line in infile:
        if( line.flag == 0):
            pos += 1
        else
            neg += 1
    print(pos)
    print(neg)
