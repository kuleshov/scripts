#!/usr/bin/env python
import argparse

# from libkuleshov.debug import keyboard
from libkuleshov.fastx import read_fasta

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--ref')
parser.add_argument('--genome')

args = parser.parse_args()

##############################################################################

ref_fasta = read_fasta(args.ref)

with open(args.genome, 'w') as genome:
	for ctg, fasta in ref_fasta.items():
		genome.write('%s\t%d\n' % (ctg, len(fasta)))
