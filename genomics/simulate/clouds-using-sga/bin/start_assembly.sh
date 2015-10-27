#!/bin/bash
set -e

#for w in {0..299}
for (( w=0; w<$1; w++ ))
do
	# echo $w
	mkdir wells/$w
	mv /home/kuleshov/metagenomica/simulate/mock-meta/simulate_clouds_by_sga/$w.shortreads.bfast.fastq wells/$w
	cd wells/$w

	SGA_PATH=/home/kuleshov/bin/asm/sga/bin

	$SGA_PATH/sga preprocess $w.shortreads.bfast.fastq | ~/lib/fastx/fastx_collapser > reads.preprocessed
	$SGA_PATH/sga index -t 16 -a ropebwt reads.preprocessed
	$SGA_PATH/sga overlap -t 16 -m 100 reads.preprocessed
	$SGA_PATH/sga assemble reads.asqg.gz

	cd -
done
