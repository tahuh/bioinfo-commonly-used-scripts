#!/usr/bin/python

"""
vcf_annotate_db.py

Annotates database id(such as dbSNP or Cosmic) to given VCF file

This program takes two input files as vcf or vcf.gz and outputs vcf file

Note. This program requires immense memory if db is large enough (not using tabix)

Author : Thomas Sunghoon Heo
"""

from __future__ import print_function
import gzip
import argparse

def fopen_gz_or_not(fname):
  handle = open(fname, "rb")
  read = handle.read(3)
  handle.seek(0)
  if read == "0x\1f\x8b\x08" :
    return gzip.GzipFile(fname, fileobj=handle)
  else:
    return handle

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--invcf", type=str, required=True, help="Input VCF. Gzip format allowed. One variant per vcf record line [no multi varianrs in one line such as GATK or MuTect gives allowed]")
parser.add_argument("-d", "--db", type=str, required=True, help="Database VCF such as dbSNP or Cosmic. Custom DB in VCF format allowed")
parser.add_argument("-o", "--outvcf", type=str, required=True, help="Output VCF file")

args = parser.parse_args()

print("[vcf_annotate_db] Start program...")
cmd = 'vcf_annodate_db.py'
cmd += " --invcf %s"%(args.invcf)
cmd += " --db %s"%(args.db)
cmd += " --outvcf %s"%(args.outvcf)

print("[vcf_annotate_db] %s"%(cmd))
print("[vcf_annotate_db] Loading database into memory. This will take time and requires huge memory if db is large")
db = {}
with fopen_gz_or_not(args.db) as dbfile:
  for line in dbfile:
    if line[0] == '#': continue
    data = line.rstrip().split("\t")
    chrom,pos,ID,ref,alt=data[0],data[1],data[2],data[3],data[4]
    for a in alt.split(','):
      K = chrom + '_' + pos + '_' + ref + '_' + a
      db[K] = ID
print ("[vcf_annotate_db] Done loading")
print ("[vcf_annotate_db] Start annotation...")
out = open(args.outvcf, "w")
nvars = 0
annots = 0
cmd_vcf_header = "##vcf_annotate_filter=%s"%(cmd)
with fopen_gz_or_not(args.invcf) as infile:
  for line in infile:
    if line[0] == '#':
      if line.startswith('#CHROM'):
        out.write(cmd_vcf_header + '\n')
        out.write(line)
      else:
        out.write(line)
      continue
    nvars += 1
    data = line.rstrip().split("\t")
    K = data[0] + '_' + data[1] + '_' + data[3] + '_' + data[4]
    try:
      db_id = db[K]
      if data[2] == '.':
        data[2] = db_id
      else:
        data[2] += ","
        data[2] += db_id
      new_line = "\t".join(data)
      annots += 1
    except KEyError:
      new_line = "\t".join(data)
    out.write(new_line)
 
 out.close()
 print("[vcf_annotate_db] %d / %d variants are annotated"%(annots,nvars))
 print("[vcf_annotate_db] Done whole process.")
