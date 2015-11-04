#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import argparse

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-o', '--out', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

bamfile = pysam.Samfile(args.bam, 'rb')
outfile = pysam.Samfile(args.out, 'wb', template=bamfile)

seen = set()
for read in bamfile:
  if read.qname not in seen:
    seen.add(read.qname)
    outfile.write(read)
