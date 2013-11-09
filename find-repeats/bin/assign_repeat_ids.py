#!/usr/bin/env python
import argparse
import os

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--psl')
parser.add_argument('--repeat_bed')

args = parser.parse_args()

##############################################################################

repeats = dict()
with os.popen('cat %s | tail -n +6' % args.psl) as psl:
	for line in psl:
		fields = line.split()
		repeat_id = fields[9]

		if repeat_id not in repeats:
			repeats[repeat_id] = list()

		t_name = fields[13]
		t_start = int(fields[15])
		t_end = int(fields[16])

		repeats[repeat_id].append((t_name, t_start, t_end))

with open(args.repeat_bed, 'w') as repeat_file:
	for id_, repeat_list in repeats.iteritems():
		if len(repeat_list) == 1: continue
		for repeat_entry in repeat_list:
			t_name, t_start, t_end = repeat_entry
			repeat_file.write('%s\t%d\t%d\t%s\n' % (t_name, t_start, t_end, id_))