#!/usr/bin/env python
import argparse

# from libkuleshov.debug import keyboard
# from libkuleshov.fastx import read_fasta

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--ctg')
parser.add_argument('--out')

args = parser.parse_args()

##############################################################################

def read_fasta(fasta_file):
	with open(fasta_file) as fasta:
		fasta_by_ctg = dict()
		num_ctgs = 0
		for line in fasta:
			if line.startswith('>'):
				num_ctgs += 1
				cur_ctg = str(num_ctgs)
				fasta_by_ctg[cur_ctg] = ""
			else:
				fasta_by_ctg[cur_ctg] += line

	return fasta_by_ctg

##############################################################################

ref_fasta = read_fasta(args.ctg)

with open(args.out, 'w') as out:
	for ctg, fasta in ref_fasta.iteritems():
		out.write('>%s\n' % ctg)
		out.write(fasta)
		
