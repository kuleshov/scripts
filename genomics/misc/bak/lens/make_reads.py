import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-p', '--variants', required=True)
parser.add_argument('-c', '--contig', required=True)
parser.add_argument('-r', '--reads', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# determine variant positions

bamfile = pysam.Samfile(args.bam, "rb")
reads = dict()
variants = dict()

n = 0
with open(args.variants) as f:
  for line in f:
    fields = line.strip().split()
    ctg, pos, alleles = fields[0], int(fields[1]), {fi.split(':')[0] for fi in fields[2:]}

    if ctg != args.contig:
      continue
    else:
      variants[pos] = alleles

# for pcol in bamfile.pileup(ctg, pos-1, pos): <-- doesn't work b/c bug in pysam
for pcol in bamfile.pileup(args.contig):
  if pcol.pos not in variants: continue
  alleles = variants[pcol.pos]

  for pread in pcol.pileups:
    if pread.indel < 0:
      # deletion
      seq = "-"
    elif pread.indel > 0:
      # we have an insertion
      seq = 'I' + pread.alignment.seq[pread.qpos+1:pread.qpos+1+pread.indel]
    else:
      # normal match
      seq = pread.alignment.seq[pread.qpos]

      q = ord(pread.alignment.qual[pread.qpos])-33
      if q < 20:
        continue

    if seq in alleles:
      read = pread.alignment.qname
      if read not in reads:
        reads[read] = dict()
      reads[read][pcol.pos] = seq

with open(args.reads, 'w') as out:
  for read in reads:
    out.write('%s\t%s' % (args.contig, read))
    for pos, seq in sorted(reads[read].iteritems()):
      if not seq.startswith('I'):
        out.write('\t%d:%s' % (pos, seq))
    out.write('\n')
