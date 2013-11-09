#!/usr/bin/env python
import argparse

# from libkuleshov.debug import keyboard
from libkuleshov.fastx import read_fasta

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--ref')
parser.add_argument('--bed')

args = parser.parse_args()

##############################################################################

ref_fasta = read_fasta(args.ref)

with open(args.bed, 'w') as bed:
	for ctg, fasta in ref_fasta.items():
		bed.write('%s\t%d\t%d\n' % (ctg, 0, len(fasta)-1))
