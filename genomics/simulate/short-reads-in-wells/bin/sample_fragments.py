#!/usr/bin/env python
import argparse
import random

from random import randint

from libkuleshov.fastx import read_fasta
from libkuleshov.misc import weighted_choice

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--ref')
parser.add_argument('--len', type=int)
parser.add_argument('--num_fragments', type=int)
parser.add_argument('--num_wells', type=int)

args = parser.parse_args()

##############################################################################

random.seed(0)

ref_fasta = read_fasta(args.ref)

for w in xrange(args.num_wells):
	with open('%d.fragments.bed' % w, 'w') as bed:
		for i in xrange(args.num_fragments):
			choices = [(ctg_name, len(ctg_fasta)) for ctg_name, ctg_fasta in ref_fasta.iteritems()]
			ctg_name = weighted_choice(choices)
			ctg_len = len(ref_fasta[ctg_name])
			
			r_start = randint(0,ctg_len-1)
			r_end = min(r_start + args.len, ctg_len-1)
			bed.write('%s\t%d\t%d\n' % (ctg_name, r_start, r_end))