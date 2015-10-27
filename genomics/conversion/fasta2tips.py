#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pysam

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', required=True)
parser.add_argument('-t', '--tips-fasta', required=True)
parser.add_argument('-l', '--tip-length', type=int, default=77)

args = parser.parse_args()

# ----------------------------------------------------------------------------

fasta = pysam.FastaFile(args.fasta)

with open(args.tips_fasta, 'w') as out:
  for ctg, length in zip(fasta.references, fasta.lengths):
    left_seq = fasta.fetch(ctg, 0, args.tip_length-1)
    right_seq = fasta.fetch(ctg, length-args.tip_length+1, length)

    out.write('>%s_LEFT\n%s\n' % (ctg, left_seq))
    out.write('>%s_RIGHT\n%s\n' % (ctg, right_seq))
