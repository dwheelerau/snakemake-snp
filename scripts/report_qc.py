#!/usr/bin/env python
import sys

# simple script to print QC results

infile = sys.argv[1]

flagin = 0
flagout = 0
with open(infile) as f:
    r1 = []
    r2 = []
    for line in f:
        bits = line.split(' ')
        for bit in bits:
            if bit.find("in1=") == 0 and len(r1) == 0:
                r1.append(bit)
                print(bit)
            elif bit.find("in2=") == 0 and len(r2) == 0:
                r2.append(bit)
                print(bit)
            elif bit.find("Input:") == 0:
                print(line)
                flagin = 1
                data = []
            elif bit.find('Result:') == 0 and flagin == 1:
                flagin = 0
                r1 = []
                r2 = []
                print("".join(data) + "\n------\n")
                data = []
            elif flagin == 1:
                data.append(bit)
            else:
                pass
