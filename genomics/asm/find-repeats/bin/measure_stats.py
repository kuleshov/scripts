#!/usr/bin/env python
import argparse

from libkuleshov.stats import n50
from libkuleshov.fastx import read_bed

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--contigs')

args = parser.parse_args()

##############################################################################

bed = read_bed(args.contigs)

lengths = [(end - start + 1) for ctg in bed.keys() for (start, end) in bed[ctg]]
print 'Contigs:', len(lengths)
print 'N50:', n50(lengths)
