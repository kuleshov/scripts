#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import itertools
import pickle
import pysam

import networkx as nx

from libkuleshov.stats import n50

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-f', '--fasta', required=True)
parser.add_argument('-l', '--rlen', type=int, default=101)

args = parser.parse_args()

# ----------------------------------------------------------------------------

bamfile = pysam.Samfile(args.bam, 'rb')
fasta = pysam.FastaFile(args.fasta)

ctg_lengths = {ctg: length for ctg, length in zip(fasta.references, fasta.lengths)}
# ctg_seq = {ctg: fasta.fetch(ctg) for ctg in fasta.references}

# ----------------------------------------------------------------------------

# compute very coarse edges: don't look at orientation, and verify number of aligned
# base pairs for only one read

# load full connections
reads = list()
edge_counts = dict()
n_ref = len(bamfile.references)
for i, ref in enumerate(bamfile.references):
  if i % 1000 == 0: print '%d/%d\n' % (i, n_ref)
  ref_len = ctg_lengths[ref]
  if ref_len > 700:
    reads = itertools.chain(
      bamfile.fetch(ref, 0,350),
      bamfile.fetch(ref, ref_len-350, ref_len)
    )
  else:
    reads = bamfile.fetch(ref)
  for read in reads:
    if not read.is_paired: continue
    if read.reference_id != read.next_reference_id:
      name1 = bamfile.getrname(read.reference_id)
      name2 = bamfile.getrname(read.next_reference_id)
      len1 = ctg_lengths[name1]
      len2 = ctg_lengths[name2]
      start1 = read.reference_start
      start2 = read.next_reference_start

      # print name1, name2, read.next_reference_start, ctg_lengths[name2]
      # print read.reference_start, read.reference_end, read.reference_end - read.reference_start + 1, len1
      # print read.cigarstring
      # print sum([l for (o,l) in read.cigartuples if o == 0])
      # print
      # assert read.reference_start <= read.reference_end

      matches = sum([l for (o,l) in read.cigartuples if o == 0])
      if matches > 100:
        if 0 <= start2 <= 270 or len2-250 <= start2 <= len2:
          edge = frozenset([name1, name2])
          if edge not in edge_counts:
            edge_counts[edge] = 0
          edge_counts[edge] += 1

pickle.dump(edge_counts, open('pe_edge_counts.pkl', 'wb'))
