set -e;
JF=/home/kuleshov/bin/jellyfish-2.1.1/bin/jellyfish
while read sampleid
do 
  cd $sampleid
  zcat *.fastq.gz > reads.tmp.fastq
  $JF bc -m 31 -s 5G -t 16 -o reads.tmp.bc reads.tmp.fastq
  $JF count -m 31 -s 5G -t 16 --bc reads.tmp.bc reads.tmp.fastq
  rm reads.tmp.fastq
  cd ..
  echo "\n"
done < $1
