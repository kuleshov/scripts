#!/usr/bin/env python

import argparse
import re

from libkuleshov.debug import keyboard

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--fasta')
parser.add_argument('--genome')
args = parser.parse_args()

##############################################################################

genome = open(args.genome, 'w')

with open(args.fasta) as f:
	fasta = f.read()

for ctg in re.finditer(r'>\w+\n[ATCG\n]+', fasta):
	m = re.search(r'^>(\w+)\n([ATCG\n]+)$', ctg.group())
	ctg_name = m.group(1)
	ctg_seq = m.group(2).replace('\n', '')

	genome.write('%s\t%d\n' % (ctg_name, len(ctg_seq)))