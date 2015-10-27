#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from libkuleshov.dna import reverse_complement

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fastg', required=True)
parser.add_argument('-a', '--asqg', required=True)
parser.add_argument('-o', '--overlap', required=True, type=int)

args = parser.parse_args()

# ----------------------------------------------------------------------------

def add_vertex(v_name, v_seq):
  prev_seq = vertices.get(v_name, None)
  if prev_seq:
    pass
    # if ori == 1:
    #   pass
    #   # this doesn't hold in general. wtf, spades!?
    #   # print len(prev_seq), len(v_seq)
    #   # print [(i, (s1, s2)) for i, (s1,s2) in enumerate(zip(prev_seq, v_seq)) if s1!=s2]
    #   # assert prev_seq == (v_seq)
    # else:
    #   pass
    #   # and this also doesn't work!!
    #   # if prev_seq != v_seq:
    #   #   print v_name
    #   #   print len(prev_seq), len(v_seq)
    #   #   print [(i, (s1, s2)) for i, (s1,s2) in enumerate(zip(prev_seq, v_seq)) if s1!=s2]
    #   # assert prev_seq == v_seq
  else:
    vertices[v_name] = v_seq

# load all the vertices and edges:
vertices = dict()
edges = list()
with open(args.fastg) as f:
  n = 0
  print 'Loading fastg'
  line1 = f.readline().strip()
  line2 = f.readline().strip()
  while line2:
    if n % 10000 == 0: print n
    if not line1.startswith('>'):
      print line1
    assert line1.startswith('>')
    line1 = line1[1:-1]
    if ':' in line1:
      # we have an edge
      origin, neighbors = line1.split(':')
      if origin.endswith("'"):
        origin_v = origin[:-1]
        origin_o = 1
      else:
        origin_v = origin
        origin_o = 0

      # add origin as vertex
      if origin_o == 0:
        v_seq = line2
      else:
        v_seq = reverse_complement(line2)
      add_vertex(origin_v, v_seq)

      for neighbor in neighbors.split(','):
        if neighbor.endswith("'"):
          edges.append( (origin_v, origin_o, neighbor[:-1], 1) )
        else:
          edges.append( (origin_v, origin_o, neighbor, 0) )
    else:
      # we have a vertex
      if line1.endswith("'"):
        v_name = line1[:-1]
        ori = 1
      else:
        v_name = line1
        ori = 0
      if ori == 0:
        v_seq = line2
      else:
        v_seq = reverse_complement(line2)

      add_vertex(v_name, v_seq)

    n += 1
    line1 = f.readline().strip()
    line2 = f.readline().strip()

# write asqg file:
placed_edges = set()
with open(args.asqg, 'w') as out:
  print 'Converting to asqg'

  print 'Converting vertices'
  for i, (v_name, v_seq) in enumerate(vertices.iteritems()):
    if i % 10000 == 0: print i
    out.write('VT\t%s\t%s\n' % (v_name, v_seq))

  print 'Converting edges'
  invalid_edges = 0
  strange_edges = 0
  for i, (v1, o1, v2, o2) in enumerate(edges):
    if i % 10000 == 0: print i
    if v1 not in vertices or v2 not in vertices:
      if v1 not in vertices: print(v1)
      if v2 not in vertices: print(v2)
      invalid_edges += 1
      continue
    if (v1,v2) in placed_edges or (v2,v1) in placed_edges:
      continue
    l1 = len(vertices[v1])
    l2 = len(vertices[v2])

    if o1 == o2:
      s1 = max(l1-args.overlap, 0)
      e1 = l1-1
      seq1 = vertices[v1]
      s2 = 0
      e2 = min(args.overlap-1, l2-1)
      seq2 = vertices[v2]
      ori = 0
    elif o1 == 0 and o2 == 1:
      s1 = max(l1-args.overlap, 0)
      e1 = l1-1
      seq1 = vertices[v1]
      s2 = max(l2-args.overlap, 0)
      e2 = l2-1
      seq2 = vertices[v2]
      ori = 1
    elif o1 == 1 and o2 == 0:
      s1 = 0
      e1 = min(args.overlap-1, l1-1)
      seq1 = vertices[v1]
      s2 = 0
      e2 = min(args.overlap-1, l2-1)
      seq2 = vertices[v2]
      ori = 1
    elif o1 == 1 and o2 == 1:
      s1 = 0
      e1 = min(args.overlap-1, l1-1)
      seq1 = vertices[v1]
      s2 = max(l2-args.overlap, 0)
      e2 = l2-1
      seq2 = vertices[v2]
      ori = 0

    if not s1 <= e1 or not s2 <= e2:
      print v1, v2
      exit()

    if ori == 0:
      # assert seq1[s1:e1+1] == seq2[s2:e2:1]
      if not seq1[s1:e1+1] == seq2[s2:e2+1]:
        print v1
        print v2
        print o1, o2
        print seq1[s1:e1+1], seq2[s2:e2+1]
        strange_edges += 1
        # continue
    else:
      # assert seq1[s1:e1+1] == reverse_complement(seq2[s2:e2:1])
      if not seq1[s1:e1+1] == reverse_complement(seq2[s2:e2+1]):
        print v1
        print v2
        print o1, o2
        print seq1[s1:e1+1], reverse_complement(seq2[s2:e2+1])
        strange_edges += 1
        # continue

    out.write('ED\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t0\n'
              % (v1, v2, s1, e1, l1, s2, e2, l2, ori))
    placed_edges.add((v1,v2))

print 'WARNING: %d edges could not be converted' % invalid_edges
print 'WARNING: %d edges correspond to funny overlaps' % strange_edges
