#!/usr/bin/env python
import argparse

# from libkuleshov.debug import keyboard
# from libkuleshov.fastx import read_fasta

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--len', type=int)
parser.add_argument('--shift', type=int)
parser.add_argument('--ref')
parser.add_argument('--fa')

args = parser.parse_args()

##############################################################################

def read_fasta(fasta_file):
	with open(fasta_file) as fasta:
		fasta_by_ctg = dict()
		num_ctgs = 0
		for line in fasta:
			if line.startswith('>'):
				num_ctgs += 1
				cur_ctg = line[1:].strip().split()[0]
				fasta_by_ctg[cur_ctg] = ""
			else:
				fasta_by_ctg[cur_ctg] += line.strip()

	return fasta_by_ctg

##############################################################################

ref_fasta = read_fasta(args.ref)

with open(args.fa, 'w') as out:
	for ctg, fasta in ref_fasta.iteritems():
		for i in xrange(0, len(fasta)-args.len, args.shift):
			out.write('>%s_%d_%d\n' % (ctg, i, i + args.len - 1))
			out.write(fasta[i:i+args.len] + '\n')
		
