import os
import sys

from libkuleshov.stats import n50

print "cat %s | tr -d '\\n'  | sed -re 's/>ctg([0-9])+/\\n/g' |     perl -nle 'print length'" % sys.argv[1]
fasta_length_pipe = os.popen("cat %s | tr -d '\\n'  | sed -re 's/>ctg([0-9])+/\\n/g' | perl -nle 'print length'" % sys.argv[1])

fasta_lengths = [int(line) for line in fasta_length_pipe]
total_len = sum(fasta_lengths)

print 'N50: %d (%d total)' % (n50(fasta_lengths), total_len)
