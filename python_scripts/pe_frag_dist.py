#!/usr/bin/python

import pysam
import argparse
from collections import defaultdict

def calculate_pe_molecule_size_exact(read1, read2):
	# Calculates actual molecule size sequenced
	# Judge if two molecules are aligned on same chromosome
	if read1.reference_name != read2.reference_name:
		return 0
	# Select which one is left most mapped read
	if read1.reference_start <= read2.reference_start:
		left = read1
		right = read2
	else:
		left = read2
		right = read1
	left_soft_or_hard_clip = 0
	right_soft_or_hard_clip = 0
	cigar_left = left.cigartuples
	if cigar_left[0][0] == 4 or cigar_left[0][0] == 5:
		left_soft_or_hard_clip += cigar_left[0][1]
	cigar_right = right.cigar_tuples
	if cigar_right[-1][0] == 4 or cigar_right[-1][0] == 5:
		right_soft_or_hard_clip += cigar_right[-1][1]
	left_start = left.reference_start
	right_start = right.reference_start
	size = right_start - left_start + 1
	qstart = read2.query_alignment_start
	qend = read2.query_alignment_end
	size += (qend - qstart + 1)
	size += left_soft_or_hard_clip
	size += right_soft_or_hard_clip
	return size
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam", type=str, required=True, help="BAM file")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output histogram data")
#parser.add_argument("-B", "--bed", type=str, required=False, help="BED file",default="")

args = parser.parse_args()

bam = pysam.AlignmentFile(args.bam)
processed = 0
pdict = {}
histo = defaultdict(int)
for record in bam.fetch():
	if record.is_unmapped :
		continue
	if record.is_supplementary:
		continue
	if record.is_secondary:
		continue
	if record.is_qcfail:
		continue
	qname = record.query_name
	try:
		pdict[qname].append(record)
	except KeyError:
		pdict[qname] = [record]
	processed += 1
	if len(pdict[qname]) == 2:
		size = calculate_pe_molecule_size_exact(pdict[qname][0],pdict[qname][1])
		pdict.pop(qname)
		histo[size] += 1
print("processed %d records\n"%(processed))
sorted_items = sorted(lambda x : x[0] , list(histo.items()))
outfile = open(args.outfile)
outfile.write("size\tcount\n")
for k, v in sorted_items:
	outfile.write(str(k) + '\t' + str(v) + '\n')
