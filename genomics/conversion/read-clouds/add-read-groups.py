#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import argparse

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-o', '--out', required=True)
parser.add_argument('-n', '--n', required=True, type=int,
                    help='number of read groups')

args = parser.parse_args()

# ----------------------------------------------------------------------------

def make_header(bamfile, n):
  """Add read group info to a header."""
  # CREATE TEMPLATE
  # Read group. Unordered multiple @RG lines are allowed.
  RG_template = { 'ID': '',           # Read group identifier. e.g., Illumina flowcell + lane name and number
                  'CN': '',           # GATK Not Required. Name of sequencing center producing the read.
                  'DS': '',           # GATK Not Required. Description
                  'DT': '',           # GATK Not Required. Date the run was produced (ISO8601 date YYYY-MM-DD or YYYYMMDD)
                  'PI': '',           # GATK Not Required. Predicted median insert size.
                  'PU': '',           # GATK Not Required. Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD).
                  'SM': '',           # Sample. Use pool name where a pool is being sequenced.
                  'PL': 'ILLUMINA'}   # Platform/technology used to produce the reads.

  new_header = bamfile.header.copy()
  new_header['RG'] = []

  # ADD INFO TO TEMPLATE
  for i in xrange(n):
    RG_template = RG_template.copy()
    RG_template['ID'] = str(i)
    new_header['RG'].append(RG_template)

  return new_header

def main():
  bamfile = pysam.Samfile(args.bam, 'rb')
  new_header = make_header(bamfile, args.n)
  # print new_header

  outfile = pysam.AlignmentFile(args.out, 'wb', header=new_header)

  for read in bamfile.fetch():
    name = read.qname
    if not name.startswith('well'):
      continue

    fields = name.split('_')
    well_id = fields[0][4:]
    new_tags = read.tags
    new_tags.append(('RG', well_id))
    read.tags = new_tags
    outfile.write(read)

  bamfile.close()
  outfile.close()

  pysam.index(args.out)

if __name__ == '__main__':
    main()
