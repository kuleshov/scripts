#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--containment', required=True)
parser.add_argument('-s', '--csv', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

# load containment
containment = dict()
with open(args.containment) as f:
  for line in f:
    fields = line.strip().split()
    ctg, chrom, start, end = fields
    if ctg not in containment:
      containment[ctg] = list()
    containment[ctg].append('%s-%s-%s' % (chrom, start, end))

out = open(args.csv, 'w')
out.write('Node name,Containment Label')
for ctg, matches in containment.iteritems():
  matches_str = ':'.join(matches)
  out.write('%s,%s\n' % (ctg, matches_str))
