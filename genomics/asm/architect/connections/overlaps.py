#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import pickle
import pysam

import nwalign as nw
from Bio.pairwise2 import align
import networkx as nx

from libkuleshov.stats import n50
from libkuleshov.dna import reverse_complement

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--containment', required=True)
parser.add_argument('-f', '--fasta', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------

fasta = pysam.FastaFile(args.fasta)
ctg_lengths = {ctg: length for ctg, length in zip(fasta.references, fasta.lengths)}
ctg_seq = {ctg: fasta.fetch(ctg) for ctg in fasta.references}

# ----------------------------------------------------------------------------

# load true intervals
ctg_ivls = dict()
with open(args.containment) as f:
  for line in f:
    fields = line.strip().split()
    if fields[1] not in ctg_ivls:
      ctg_ivls[fields[1]] = list()
    ctg_ivls[fields[1]].append((fields[2], int(fields[3]), int(fields[4])))

if os.path.isfile('./connected_components.pkl'):
  true_overlaps = pickle.load(open('./true_overlaps.pkl', 'rb'))
  cc_list = pickle.load(open('./connected_components.pkl', 'rb'))
else:
  print 'computing overlaps'

  # compute true overlaps between contigs
  true_overlaps = set()
  for i, (ctg1, ivls1) in enumerate(ctg_ivls.iteritems()):
    if i % 100 == 0: print '%d/%d' % (i, len(ctg_ivls))
    for ctg2, ivls2 in ctg_ivls.iteritems():
      if any(ivl_dist(i1,i2) < 200 for i1 in ivls1 for i2 in ivls2):
        true_overlaps.add(frozenset([ctg1, ctg2]))

  pickle.dump(true_overlaps, open('true_overlaps.pkl', 'wb'))

  print 'computing components'

  # compute connected components
  G = nx.Graph()
  G.add_nodes_from(fasta.references)
  for i, ovl1 in enumerate(true_overlaps):
    if i % 1000 == 0: print '%d/%d' % (i, len(true_overlaps))
    ol = list(ovl1)
    if len(ol) > 2:
      assert False
    elif len(ol) < 2:
      continue
    G.add_edge(ol[0],ol[1])
  connected_components = nx.connected_components(G)

  cc_list = [cc for cc in connected_components]
  pickle.dump(cc_list, open('connected_components.pkl', 'wb'))

component_lengths = [sum([ctg_lengths[ctg] for ctg in component]) 
                     for component in cc_list]

print sum(fasta.lengths)
print sum(component_lengths)
print n50(component_lengths)

# ----------------------------------------------------------------------------
# compute overlaps

def myalign(s1,s2):
  return align.globalms(s1,s2,2,-1,-2,-0.2,penalize_end_gaps=False)[0]

def ctgs_overlap(ctg1, ctg2):
  head_seq1 = ctg_seq[ctg1][:100]
  tail_seq1 = ctg_seq[ctg1][-100:]
  head_seq2 = ctg_seq[ctg2][:100]
  tail_seq2 = ctg_seq[ctg2][-100:]

  hh_align = myalign(head_seq1, reverse_complement(head_seq2))
  hh_score = hh_align[2]
  ht_align = myalign(head_seq1, tail_seq2)
  ht_score = ht_align[2]
  th_align = myalign(tail_seq1, head_seq2)
  th_score = th_align[2]
  tt_align = myalign(tail_seq1, reverse_complement(tail_seq2))
  tt_score = tt_align[2]

  # hh_align = nw.global_align(head_seq1, reverse_complement(head_seq2))
  # hh_score = len([x for x,y in zip(*hh_align) if x == y])
  # ht_align = nw.global_align(head_seq1, tail_seq2)
  # ht_score = len([x for x,y in zip(*ht_align) if x == y])
  # th_align = nw.global_align(tail_seq1, head_seq2)
  # th_score = len([x for x,y in zip(*th_align) if x == y])
  # tt_align = nw.global_align(tail_seq1, reverse_complement(tail_seq2))
  # tt_score = len([x for x,y in zip(*tt_align) if x == y])

  scores = (hh_score, ht_score, th_score, tt_score)
  aligns  =  (hh_align, ht_align, th_align, tt_align)

  best_i = max([0,1,2,3],key=lambda x: scores[x])
  # print best_i
  # print scores[best_i]
  # print aligns[best_i][0]
  # print aligns[best_i][1]

  if scores[best_i] > 100:
    return True
    # overlaps.add(frozenset([ctg1,ctg2]))

# overlaps = set()
# ctgs = list(ctg_seq.keys())
# for i, ctg1 in enumerate(ctgs):
#   # if i % 20 == 0: print '%d/%d' % (i, len(ctgs))
#   print '%d/%d' % (i, len(ctgs))
#   for ctg2 in ctgs[i:]:
#     if ctg1 == ctg2: continue
#     if ctgs_overlap(ctg1, ctg2):
#       overlaps.add(frozenset([ctg1,ctg2]))

# ----------------------------------------------------------------------------
# verify true overlaps:

validated_overlaps = set()
for i, ovl in enumerate(true_overlaps):
  if i % 1000 == 0: print i
  if len(ovl) == 1:
    validated_overlaps.add(ovl)
    continue
  ctg1, ctg2 = list(ovl)
  # if ctg_lengths[ctg1] < 1000 or ctg_lengths[ctg2] < 1000:
  #   continue
  # print
  # print ctg1, ctg_lengths[ctg1], ctg_ivls[ctg1]
  # print ctg2, ctg_lengths[ctg2], ctg_ivls[ctg2]
  if ctgs_overlap(ctg1, ctg2):
    validated_overlaps.add(ovl)


print len(true_overlaps), len(validated_overlaps)

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
