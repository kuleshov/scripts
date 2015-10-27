def read_fasta(ctg_file, parse=lambda x: '_'.join(x[1:].strip().split())):
	ctgs = dict()
	with open(ctg_file) as ctg:
		line = ctg.readline()
		while line:
			if line.startswith('>'):
				ctg_name = parse(line)
				ctgs[ctg_name] = ''
			else:
				ctgs[ctg_name] += line.strip()
			line = ctg.readline()

	return ctgs

def read_bed(bed_path):
	with open(bed_path) as bed_file:
		bed = dict()
		for line in bed_file:
			ctg, start, end = line.strip().split()[:4]
			if ctg not in bed:
				bed[ctg] = list()
			bed[ctg].append((int(start),int(end)))

	return bed

def to_80_chars(line):
	num_lines = len(line) / 80
	out = ''
	for i in xrange(num_lines):
		out += (line[i*80:(i+1)*80] + '\n')

	out += line[num_lines*80:] + '\n'

	return out
