import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-c', '--contig', required=True)
parser.add_argument('-p', '--pos', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# determine variant positions

bamfile = pysam.Samfile(args.bam, "rb")
posfile = open(args.pos, "w")

n = 0
for pcol in bamfile.pileup(args.contig):
  alleles = dict()
  for pread in pcol.pileups:
    if pread.alignment.query_qualities[pread.query_position] < 15: continue
    if pread.indel != 0: continue
    if pread.is_del: continue
    if pread.indel < 0 or pread.is_del:
      # deletion
      seq = "-"
    elif pread.indel > 0:
      # we have an insertion
      seq = 'I' + pread.alignment.seq[pread.query_position:pread.query_position+1+pread.indel]
    else:
      # normal match
      seq = pread.alignment.query_sequence[pread.query_position]

    if seq not in ('N',):
      if seq not in alleles:
        alleles[seq] = 0
      alleles[seq] += 1

  total_alleles = sum(alleles.values())
  threshold = max(3, 0.1*total_alleles)
  significant_variants = {k:v for (k,v) in alleles.iteritems() if v >= threshold}
  if len(significant_variants) >= 2:
    posfile.write('%s\t%d' % (args.contig, pcol.pos))
    for b, n in sorted(significant_variants.iteritems(), key=lambda x: x[1], reverse=True):
      posfile.write('\t%s:%d' % (b, n))
    posfile.write('\n')

  n += 1
  if n % 100000 == 0:
    print n



