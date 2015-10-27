#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--inp', required=True)
parser.add_argument('-o', '--out', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# permute bases

outfile = open(args.out, 'w')
n = 1

with open(args.inp) as infile:
  line1 = infile.readline()
  line2 = infile.readline()

  while line2:
    line2a = list(line2)
    idx = [i for (i,c) in enumerate(line2a) if c == 'N']
    for i in idx:
      line2a[i] = random.choice(('A', 'T', 'C', 'G'))

    outfile.write(line1)
    outfile.write(''.join(line2a))

    line1 = infile.readline()
    line2 = infile.readline()

    n += 1
    if (n % 100000) == 0: print n
