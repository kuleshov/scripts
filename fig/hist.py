#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from matplotlib import pyplot as plt

# ----------------------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument('-o', '--out', required=True)
parser.add_argument('-x', '--max', type=int)
parser.add_argument('-n', '--min', type=int)
parser.add_argument('-t', '--title')

args = parser.parse_args()

# ----------------------------------------------------------------------------

X = [float(line) for line in sys.stdin]
plt.hist(X,bins=50)
# plt.show()
max_val, min_val = max(X), min(X)
if args.max:
  max_val = min(max(X),args.max)
if args.min:
  min_val = min(min(X), args.min)
if args.title:
  plt.title(args.title)
plt.xlim(min_val, max_val)
plt.savefig(args.out)
