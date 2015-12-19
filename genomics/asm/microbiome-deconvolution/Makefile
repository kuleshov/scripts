# RECIPE FOR CREATING CONTIG CLUSTERS

#------------------------------------------------------------------------------
# config

PROCESSES=10

# execution directories
SCRIPTDIR:=/home/kuleshov/work/mf-asm/src
export ROOT:=/home/kuleshov/brian-yu

# program paths
PARALLEL:=parallel
export SPADES:=/home/kuleshov/bin/asm/SPAdes-3.5.0-Linux-compiled/bin/spades.py
export ROOT:=/home/kuleshov/brian-yu
export BWA:=/home/kuleshov/bin/bwa
export SAMTOOLS:=/sw/bin/samtools

#------------------------------------------------------------------------------
# recipe

all:
	make asm.done
	make align.done
	make cluster.done

asm.done:
	# copy names of barcodes we will consider (might have to regenerate)
	cp $(SCRIPTDIR)/asm/barcodes.txt .
	# assemble each barcode using spades
	bash $(SCRIPTDIR)/asm/asm.sh | $(PARALLEL) -j $(PROCESSES)
	touch $@

align.done:
	# # index contigs for bwa
	bash $(SCRIPTDIR)/align/index.sh | $(PARALLEL) -j $(PROCESSES)
	# # align reads from every well to every contig file
	bash $(SCRIPTDIR)/align/align.sh | $(PARALLEL) -j $(PROCESSES)
	# merge files from the same barcode (and index the result)
	bash $(SCRIPTDIR)/align/merge.sh | $(PARALLEL) -j $(PROCESSES)
	touch $@

cluster.done:
	# generate profiles
	find . | grep bam | grep -v bai | \
	  parallel -j $(PROCESSES) python $(SCRIPTDIR)/cluster/create_profiles.py \
	    -b {} \
	    -p {.}.profile \
	    -c 1 100
	# merge profiles
	python $(SCRIPTDIR)/cluster/merge.py -b barcodes.txt -d . -o merged.profile
	# take only the 10kbp reads
	python $(SCRIPTDIR)/cluster/filter.py -p merged.profile -l 10000 > merged.10k.profile
	# cluster profiles
	python $(SCRIPTDIR)/cluster/cluster.py -p merged.10k.profile -c clusters.txt -s singletons.txt -t 10
	# touch $@
