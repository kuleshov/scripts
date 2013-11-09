#!/usr/bin/env python
import argparse
import sys
import random

from random import randint
from libkuleshov.fastx import read_fasta
from libkuleshov.debug import keyboard

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--len', type=int)
parser.add_argument('--cov', type=float)
parser.add_argument('--even', action='store_true')
parser.add_argument('--meta')
parser.add_argument('--ref')
parser.add_argument('--qual')
parser.add_argument('--fasta')

args = parser.parse_args()

##############################################################################

ctg_by_org = dict()
weight_by_org = dict()
with open(args.meta) as meta:
	for line in meta:
		fields = line.split()
		ctg = int(fields[0])
		organism = fields[-1]
		weight = fields[4]

		if organism not in ctg_by_org:
			ctg_by_org[organism] = [ctg]
			weight_by_org[organism] = weight
		else:
			ctg_by_org[organism].append(ctg)

##############################################################################

def weighted_choice(choices):
   total = sum(w for c, w in choices)
   r = random.uniform(0, total)
   upto = 0
   for c, w in choices:
      if upto + w >= r:
         return c
      upto += w
   print choices
   assert False, "Shouldn't get here"

##############################################################################

random.seed(0)

parse_fn = lambda x: int(x[1:].strip().split('|')[1])
ctg_fasta = read_fasta(args.ref, parse=parse_fn)

genome_lengths = [sum([len(ctg_fasta[ctg]) for ctg in org_contigs]) for org_contigs in ctg_by_org.values()]
organisms = ctg_by_org.keys()
total_sequence = sum(genome_lengths)

num_reads = int(total_sequence * args.cov / args.len)
num_organisms = len(ctg_by_org.keys())

positions_covered = dict()
for o in organisms:
	positions_covered[o] = dict()
	for ctg in ctg_by_org[o]:
		positions_covered[o][ctg] = set()

fasta = open(args.fasta, 'w')
qual = open(args.qual, 'w')

total_read_len = 0

for i in xrange(num_reads):
	o = random.choice(organisms)
	choices = [(ctg, len(ctg_fasta[ctg])) for ctg in ctg_by_org[o]]
	ctg = weighted_choice(choices)

	r = randint(0,len(ctg_fasta[ctg])-1)
	read = ctg_fasta[ctg][r:r+args.len]
	total_read_len += len(read)
	positions_covered[o][ctg].update(range(r,r+len(read)))
	fasta.write('>simread_%d_%d\n' % (i, len(read)))
	fasta.write(read + '\n')

	qual.write('>simread_%d_%d\n' % (i, len(read)))
	qual.write(' '.join(['99'] * len(read)) + '\n')

print 'EFFECTIVE COVERAGE: %f (%d / %d)' % (float(total_read_len) / float(total_sequence), total_read_len, total_sequence) 

##############################################################################
## COMPUTE HOLES

total_contigs = 0
for o in organisms:
	for ctg in ctg_by_org[o]:
		uncovered = set(range(len(ctg_fasta[ctg]))) - positions_covered[o][ctg]
		uncovered.update((-1, len(ctg_fasta[ctg])))
		uncovered_list = sorted(list(uncovered))
		contig_starts = [u for i, u in enumerate(uncovered_list[1:]) if (uncovered_list[i+1] - uncovered_list[i]) > 1]
		num_contigs = len(contig_starts)
		print 'ORG: %s CTG: %s CONTIGS: %d' % (o, ctg, num_contigs)
		total_contigs += num_contigs
print 'TOTAL: %d contigs' % total_contigs

