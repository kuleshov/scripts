#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pickle

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--scaffolds', required=True)
parser.add_argument('-e', '--edges', required=True)
parser.add_argument('-t', '--tsv', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

scaff_map = dict()
with open(args.scaffolds) as f:
  for line in f:
    s,c = line.strip().split()
    scaff_map[s] = c

tsv = open(args.tsv, 'w')
edge_counts = dict()

with open(args.edges) as f:
  for line in f:
    s1, s2, conn, cnt = line.strip().split()
    # if s1 not in scaff_map or s2 not in scaff_map:
    #   continue
    c1, c2 = scaff_map[s1], scaff_map[s2]
    conn1, conn2 = conn
    if conn1 == conn2:
      strand = 'S'
    else:
      strand = 'R'
    tsv.write('%s\t%s\t%s\t%s\t%s\t%s\t100\n'
        % (c1, c2, conn1, conn2, strand, cnt))
    edge_counts[frozenset([c1,c2])] = int(cnt)

pickle.dump(edge_counts, open('besst_edge_counts.pkl', 'wb'))


