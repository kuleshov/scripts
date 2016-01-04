#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pysam

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bed', required=True)
parser.add_argument('-c', '--containment', required=True)
parser.add_argument('-r', '--ref-fasta', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

fasta = pysam.FastaFile(args.ref_fasta)
ref_names = [name.split('_')[0] for name in fasta.references]
ref_map = {name : i for i,name in sorted(enumerate(ref_names),
                               key=lambda x: fasta.lengths[x[0]],
                               reverse=True) }
def ivl_dist(ivl1, ivl2):
  chr1, chr2 = ivl1[0], ivl2[0]
  i1, i2 = ivl1[1:], ivl2[1:]
  if chr1 != chr2:
    return float('inf')
  else:
    if i1[0] <= i2[0] <= i1[1]:
      return 0
    elif i1[0] <= i2[1] <= i1[1]:
      return 0
    else:
      return min(abs(i2[1] - i1[0]), abs(i2[0] - i1[1]))

intervals = dict()
with open(args.bed) as f:
  for line in f:
    fields = line.strip().split()
    ctg = fields[3]
    if ctg not in intervals: intervals[ctg] = list()
    # ivl_str = '%s-%s-%s' % (fields[0], fields[1], fields[2])
    ivl = fields[:3]
    intervals[ctg].append((fields[0], int(fields[1]), int(fields[2])))

with open(args.containment, 'w') as out:
  for ctg, ivls in intervals.iteritems():
    # first determine major intervals:
    # major_ivls = [i for i in ivls if i[2]-i[1] > 250]
    # # take the intervals that are within 200bp of a major interval
    # valid_ivls = [i for i in ivls if any(ivl_dist(i,i_major) < 200 for i_major in major_ivls)]
    # sort intervals by chromosome:
    chroms = set([i[0] for i in ivls])
    # choose main chromosome:
    chrom = sorted(list(chroms), key=lambda x: sum([abs(i[2]-i[1]) for i in ivls if i[0] == x]),
                  reverse=True)[0]
    valid_ivls = sorted([i for i in ivls if i[0] == chrom], key=lambda x:x[1])
    grouped_ivls = list()
    curr_ivl = list(valid_ivls[0])
    for ivl in valid_ivls[1:]:
      if ivl[2] - curr_ivl[2] < 750:
        curr_ivl[2] = ivl[2]
      else:
        grouped_ivls.append(curr_ivl)
        curr_ivl = list(ivl)
    grouped_ivls.append(curr_ivl)

    for ivl in grouped_ivls:
      out.write('R\t%s\t%d\t%d\t%d\n' % (ctg, ref_map[ivl[0]], ivl[1], ivl[2]))
