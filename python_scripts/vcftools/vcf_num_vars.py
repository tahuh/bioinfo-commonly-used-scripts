#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", type=str, required=True, help="VCF file")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file")
parser.add_argument("-s", "--sample",default="sample", type=str, help="Sample name")
args = parser.parse_args()

with open(args.outfile, "w") as outfile:
	n = 0
	with open(args.vcf) as vcf:
		for line in vcf:
			if line[0] == '#': continue
			n += 1
		outfile.write("sample\tnum_vars\n")
		outfile.write(args.sample + '\t' + str(n) + '\n')

