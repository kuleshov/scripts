#!/usr/bin/env bash

# make well-barcode map
ls | grep LR |grep fastq | cut -d '_' -f 3 | sort | uniq | nl  > wells.map


# concatenate all the reads
#echo "" > reads.gz

while read line
do
  a=( $line )
  # write wells for R1 reads:
  zcat LR6000033-DNA_A01-LRAAD-01_${a[1]}_L006_R1_001.fastq.gz \
    | awk -v W="${a[0]}" '{if ($0 ~ /^@/){print "well"W"_"$1"/1"}else{print $0}}' \
    | awk 'NR % 16 <= 3' \
    | pigz \
    >> reads.1.gz

  # write wells for R2 reads:
  zcat LR6000033-DNA_A01-LRAAD-01_${a[1]}_L006_R2_001.fastq.gz \
  | awk -v W="${a[0]}" '{if ($0 ~ /^@/){print "well"W"_"$1"/2"}else{print $0}}' \
  | awk 'NR % 16 <= 3' \
  | pigz \
  >> reads.2.gz
done < wells.map
