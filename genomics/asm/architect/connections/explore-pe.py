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

edge_counts = pickle.load(open('./besst_edge_counts.pkl', 'rb'))
# edge_counts = pickle.load(open('./pe_accurate_edge_counts.pkl', 'rb'))
true_overlaps = pickle.load(open('./true_overlaps.pkl', 'rb'))

# true_ovl_list = [list(o) for o in true_overlaps]
# true_overlaps = {frozenset(o) for o in true_ovl_list
#                 if len(o) == 2
#                 and ctg_lengths[o[0]] > 700 and ctg_lengths[o[1]] > 700}

overlaps = {e for e,c in edge_counts.iteritems() if c >= 3}

validated_overlaps = overlaps & true_overlaps

print len(overlaps), len(validated_overlaps), len(true_overlaps)

# validated_overlaps = true_overlaps

# compute connected components
G = nx.Graph()
G.add_nodes_from(fasta.references)
for i, ovl1 in enumerate(validated_overlaps):
  if i % 1000 == 0: print '%d/%d' % (i, len(validated_overlaps))
  ol = list(ovl1)
  if len(ol) > 2:
    assert False
  elif len(ol) < 2:
    continue
  G.add_edge(ol[0],ol[1])
connected_components = nx.connected_components(G)

cc_list = [cc for cc in connected_components]

component_lengths = [sum([ctg_lengths[ctg] for ctg in component])
                     for component in cc_list]

print sum(component_lengths)
print n50(component_lengths)
