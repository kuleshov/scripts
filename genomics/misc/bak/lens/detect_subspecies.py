import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bam', required=True)
parser.add_argument('-c', '--contig', required=True)
parser.add_argument('-r', '--reads', required=True)
parser.add_argument('-k', '--clusters', required=True)

args = parser.parse_args()

# ----------------------------------------------------------------------------
# constants

OVERLAP_CUTOFF = 2            # min bases shared between two merged reads
SIMILARITY_CUTOFF = 1.0       # percentage of the bases of one read that
                              # must be the same across both reads
COV_CUTOFF = 3               # min avergae coverage of a block
COV_CUTOFF2 = (2, 0.75)       # at least y of block must have coverage x

# ----------------------------------------------------------------------------
# load reads

reads = dict()
with open(args.reads) as f:
  for line in f:
    fields = line.strip().split()
    ctg, read = fields[0], fields[1]
    if ctg != args.contig:
      continue

    alleles = {int(fi.split(':')[0]):fi.split(':')[1] for fi in fields[2:]}
    reads[read] = alleles

# ----------------------------------------------------------------------------
# classes/functions used in main algorithm

class Profile:
  def __init__(self, profile):
    # _profile is a dict of the form {pos:allele}
    assert profile
    self._profile = profile

  def common_pos(self, other):
    return [allele for allele in self._profile if allele in other._profile]

  def similar(self, profile):
    n_equal = len([x for x in self if x in profile])
    n_common = len(self.common_pos(profile))
    if n_common >= OVERLAP_CUTOFF and (n_equal >= SIMILARITY_CUTOFF*n_common):
      return True
    else:
      return False

  def __len__(self):
    return len(self._profile)

  def __iter__(self):
    return self._profile.iteritems()

  def __contains__(self, x):
    return (self._profile.get(x[0], None) == x[1])

  def __eq__(self, other):
    fz0 = frozenset(self._profile.items())
    fz1 = frozenset(other._profile.items())
    return fz0 == fz1

  def __hash__(self):
    return frozenset(self._profile.items()).__hash__()

  def update(self, profile):
    new_profile = dict(self._profile.items() + profile._profile.items())
    return Profile(new_profile)

  def start(self):
    return sorted(self._profile.keys())[0]

  def end(self):
    return sorted(self._profile.keys(), reverse=True)[0]

  def __str__(self):
    sorted_items = sorted((x for x in self._profile.iteritems()), key=lambda x: x[0])
    return ''.join(x[1] for x in sorted_items)

  def string_over(self, positions):
    return ''.join([self._profile.get(pos, '?') for pos in sorted(positions)])

  def delete(self, pos):
    del self._profile[pos]

def write_clusters(clusters, n):
  cluster_file = open(args.clusters, 'a')
  sorted_items = sorted(clusters.iteritems(), key=lambda x: (x[0].start(), x[0].end()))
  region_positions = {p for profile, R in sorted_items for (p,x) in profile}
  for profile, R in sorted_items:
    if len(R) < 3:
      continue

    profile_str = profile.string_over(region_positions)
    #print profile_str

    block_positions = sorted(profile._profile.keys())
    coverages = {pos:0 for pos in block_positions}
    for read in R:
      for pos in block_positions:
        if pos in reads[read]: coverages[pos] += 1
    coverage = float(sum(coverages.values())) / len(coverages)

    for read in R:
      read_seq = [reads[read].get(pos, '?') for pos in sorted(profile._profile.keys())]
      #print ''.join(read_seq)

    # print coverages if it has at least two coverage over 75% of its length:
    if coverage < COV_CUTOFF or (len([c for c in coverages.values() if c >= COV_CUTOFF2[0]]) < COV_CUTOFF2[1]*len(block_positions)):
      #print 'REJECTED\n'
      continue
    else:
      #print
      pass

    for pos in block_positions:
      if coverages[pos] < 2:
      	profile.delete(pos)

    cluster_name = args.contig + '-' + str(n)
    cluster_file.write('%s\t%d\t%d\t%s\t%f\t%s\t' % (args.contig, profile.start(), profile.end(), cluster_name, coverage, profile_str))
    for p, a in profile:
      cluster_file.write('%d:%s,' % (p,a))
    for read in R:
      cluster_file.write('\t%s' % read)
      read_seq = [reads[read].get(pos, '?') for pos in sorted(profile._profile.keys())]

    #print
    cluster_file.write('\n')

  cluster_file.close()
   
# ----------------------------------------------------------------------------
# determine groups of reads

bamfile = pysam.Samfile(args.bam, 'rb')

# touch cluster file
cluster_file = open(args.clusters, 'w')
cluster_file.close()

curr_reads = set()
clusters = dict()

poslist = {p for r in reads for p in reads[r]}
n = 0
# move from left to right and add reads to current group
for pos in sorted(poslist):
  #print pos
  # fetch new reads at this position:
  reads_at_pos = {r.qname for r in bamfile.fetch(args.contig, pos, pos+1)}
  old_reads = reads_at_pos & curr_reads
  new_reads = reads_at_pos - curr_reads

  if curr_reads and not old_reads:
    write_clusters(clusters,n)
    n += 1
    clusters = dict()
    curr_reads = set()

  # assign new reads to clusters, or form new clusters
  for r in new_reads:
    if r in reads:
      read_profile = Profile(reads[r])
    else:
      continue
    for cluster_profile in clusters:
      if read_profile.similar(cluster_profile):
        #print read_profile, ' similar to ', cluster_profile
        new_profile = cluster_profile.update(read_profile)
        clusters[new_profile] = clusters.pop(cluster_profile)
        clusters[new_profile].add(r)
        break
    else:
      clusters[read_profile] = set([r])
      #print 'new block', read_profile

  curr_reads.update(new_reads)

# write last clusters:
write_clusters(clusters, n)
