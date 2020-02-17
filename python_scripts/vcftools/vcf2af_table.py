#!/usr/bin/python
"""
vcf2af_table.py

Calculates AF using DP and AD field of vcf file given

Generates TSV file with 9 columns
SampleID,chrom,pos,id,ref,alt,totalDP,altDP,AF

Author : Thomas Sunghoon Heo
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-v" , "--vcf", type=str, required=True, help="VCF file")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file")
parser.add_argument("--af-as-percent", action="store_true", default=False, help="Calculate as Percentage",dest="as_percent")
args = parser.parse_args()
outfile = open(args.outfile, "w")
outfile.write("SampleID\tchrom\tpos\tid\tref\talt\ttotalDP\taltDP\tAF\n")

idx2samp = {}
with open(args.vcf) as vcf:
	for line in vcf:
		if line[0] == '#':
			if line.startswith("#CHROM"):
				data = line.rstrip().split("\t")
				samples = data[9:]
				for i, samp in enumerate(samples):
					idx2samp[i] = samp
			continue
		data = line.rstrip().split("\t")
		chrom,pos,id,ref,alts=data[0],data[1],data[2],data[3],data[4],data[5].split(",")
		tmp = {}
		fmt_str = data[8].split(":")
		fmt_values = list(map(lambda x : x.split(":"), data[9:]))
		for i, zstr in enumerate(fmt_str):
			tmp[zstr] = []
			for v in fmt_values:
				value = v[i]
				tmp[zstr].append(v)
		for i in range(len(fmt_values)):
			samp = idx2samp[i]
			totalDP = tmp['DP'][i] # i-th sample's total depth
			ADs = tmp['AD'][i].split(",") # i-th sample's allele depth
			for j, a in enumerate(alts):
				ad = int(ADs[j+1])
				af = float(ad) / float(totalDP)
				if args.as_percent:
					af = 100 * af
				outfile.write(samp + '\t' + chrom + '\t' + pos + '\t' + id + '\t' + ref + '\t' + a + '\t' + totalDP + '\t' + str(ad) + '\t' + str(af) + '\n')
outfile.close()
