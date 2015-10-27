#!/usr/bin/env python
import argparse
import os

from libkuleshov.stats import n50

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--contigs')

args = parser.parse_args()

##############################################################################

with os.popen("cat %s | grep '>'" % args.contigs) as contigs:
	lengths = list()
	for line in contigs:
		name, length, unk = line[1:].strip().split()
		lengths.append(int(length))

print 'N50:', n50(lengths)