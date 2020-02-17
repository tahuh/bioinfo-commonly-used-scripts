#!/usr/bin/python

from __future__ import print_function
import sys

bedname = sys.argv[1]
size = 0
with open(bedname) as bedfile:
	for line in bedfile:
		data = line.rstrip().split("\t")
		chrom,start,end=data[0],int(data[1]),int(data[2])
		size += (end-start+1)
print(size)
