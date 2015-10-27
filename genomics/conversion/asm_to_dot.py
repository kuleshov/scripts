#!/usr/bin/env python

import argparse
import re

from libkuleshov.debug import keyboard

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--asm')
parser.add_argument('--connections')
parser.add_argument('--nodes')
args = parser.parse_args()

##############################################################################

# load all contig and link messages

messages = list()
with open(args.asm) as asm:
	braces_stack = 0
	for line in asm:
		if line.startswith('{'):
			braces_stack += 1
			if braces_stack == 1:
				message_type = line[1:4]
				message_contents = ""
		elif line.startswith('}'):
			assert braces_stack > 0
			if braces_stack == 1:
				if message_type == 'CCO' or message_type == 'CLK':
					assert len(message_contents) > 0
					messages.append((message_type, message_contents))
			braces_stack -= 1
		else:
			if message_type == 'CCO' or message_type == 'CLK':
				message_contents += line


# extract contig and link information
contigs = dict()
links = list()

# compile regexps
contig_id_re = re.compile(r'acc:\((.+),.+\)\n')
contig_len_re = re.compile(r'len:(\d+)\n')
contig_seq_re = re.compile(r'cns:\n([ATCG\-\n]+)\n.\n')
contig_reads_re = re.compile(r'npc:(\d+)\n')
contig_complexity_re = re.compile(r'(..+)\1{9,}')

link_ids_re = re.compile(r'co1:(.+)\nco2:(.+)\n')
link_num_re = re.compile(r'num:(\d+)\n')
link_num_re = re.compile(r'mea:(-?\d+\.\d+)\n')

for message_type, message_contents in messages:
	if message_type == 'CCO':
		contig_id = contig_id_re.search(message_contents).group(1)
		# can also do re.search(r'acc:\(.+,(.+)\)\n', contents).group(1)
		contig_len = int(contig_len_re.search(message_contents).group(1))
		contig_seq = contig_seq_re.search(message_contents).group(1).replace('\n', '')

		contig_complexity = 0
		m = contig_complexity_re.search(contig_seq)
		if m:
			if float(len(m.group())) / len(contig_seq) > 0.15:
				contig_complexity = 1
				print contig_seq
				print m.group(), m.span()
		
		contigs[contig_id] = dict(len=contig_len, seq=contig_seq, compl=contig_complexity)
	elif message_type == 'CLK':
		m = link_ids_re.search(message_contents)
		co1, co2 = m.group(1), m.group(2)
		links.append((co1, co2))

connections = open(args.connections, 'w')
for co1, co2 in links:
	connections.write('%s\t0\t%s\n' % (co1, co2))

nodes = open(args.nodes, 'w')
for co in contigs:
	nodes.write('%s\t%d\t%s\t%d\n' % (co, contigs[co]['len'], contigs[co]['seq'], contigs[co]['compl']))