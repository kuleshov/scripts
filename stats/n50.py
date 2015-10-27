#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from libkuleshov.stats import n50
import sys

def n50(l,total_length=None):
    if not total_length: total_length = sum(l)
    half_length = total_length / 2
    length = 0

    # l.sort(reverse=True)
    l.sort()
    # print l

    while length < half_length and l:
        x = l.pop()
        length += x
        print "%d/%d" % (length, total_length), x

    return x

L=[int(x) for x in sys.stdin]
if len(sys.argv) > 1:
  print n50(L, float(sys.argv[1]))
else:
  print n50(L)
