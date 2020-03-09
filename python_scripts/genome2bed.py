#!/usr/bin/python

from __future__ import print_function
import argparse
import gzip

def format_checker(fqname) :
	F = open(fqname,"rb")
	bytes = F.read(3)
	if bytes == "\x1f\x8b\x08" :
		F.seek(0)
		return gzip.GzipFile(fqname , fileobj=F)
	else:
		F.seek(0)
		return F

class FAParser(object):
	def __init__(self, faname=None):
		if not faname :
			self.faname = None
		else:
			self.faname = faname

	def open(self, faname=None) :
		if self.faname == None:
			if faname == None:
				raise AttributeError("Fasta File name must be specified")

			else:
				self.faname = faname

		self.fafile = open(self.faname)
	def rewind(self):
		self.fafile.seek(0)
	def close(self) :
		self.fafile.close()

	def parse(self) :
		id = ''
		desc = ''
		seq = ''
		seq_trail = []
		for line in self.fafile:
			if line.startswith(">") :
				if seq_trail :
					yield id , desc , "".join(seq_trail)
				id = line.rstrip().split()[0][1:]
				desc = line.rstrip()[1:]
				seq_trail = []

			else:
				seq_trail.append(line.rstrip())
		if seq_trail:
			yield id , desc , "".join(seq_trail)
			
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genome", type=str, required=True, help="Genome FASTA file(accepts gzip or bgzip files)")
parser.add_argument("-o", "--output", type=str, required=True, help="Output BED file(zero based half open coordinate)")

args = parser.parse_args()
bed = open(args.output, "w")
print("Running task")

fasta = FAParser(args.genome)
fasta.open()
for name, desc, seq in fasta.parse():
	start = 0
	end = len(seq)
	bed.write(name + '\t' + str(start) + '\t' + str(end) + '\n')

print("Done task")
