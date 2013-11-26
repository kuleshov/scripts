## ARCHITECT PIPELINE

# binaries
ARCH:=~/architect/architect.py
SGA:=~/bin/asm/sga/bin/sga
FASTX:=~/lib/fastx

# input data
READS:=/home/kuleshov/metagenomica/simulate/rs/simulate_short_reads/shortreads.fastq

# parameters
CHECKPOINTS:=./checkpoints

##############################################################################

clean:
	-rm reads.*
	-rm default*

run:
	[ -e $(CHECKPOINTS) ] || mkdir -p $(CHECKPOINTS)
	make 3.resolve_repeats.done

##############################################################################
## SGA
# This run the BWT and overlap components of SGA.
# Currently, it also starts the assembly and takes the overlap graph as it
# is after removing contained vertices

reads.preprocessed: $(READS)
	$(SGA) preprocess $< | $(FASTX)/fastx_collapser > $@.tmp
	mv $@.tmp $@

reads.rbwt: reads.preprocessed
	$(SGA) index -t 16 -a ropebwt $<

reads.asqg.gz: reads.preprocessed reads.rbwt
	$(SGA) overlap -t 16 -m 100 $<

reads-contigs.fa: reads.asqg.gz
	$(SGA) assemble -o reads $<

reads.pre_simplify: reads-contigs.fa

##############################################################################
## ARCHITECT

0.load.done: reads.pre_simplify
	$(ARCH) load \
		--asqg $< \
		--out $(CHECKPOINTS)/0.loaded
	touch $@

1.contract.done: 0.load.done
	$(ARCH) contract \
		--in $(CHECKPOINTS)/0.loaded \
		--out $(CHECKPOINTS)/1.contracted
	touch $@

2.prune.done: 1.contract.done
	$(ARCH) prune \
		--in $(CHECKPOINTS)/1.contracted \
		--out $(CHECKPOINTS)/2.pruned \
		> 2.prune.log
	touch $@

3.resolve_repeats.done 3.resolve.log 3.resolved.stats 3.resolved.dot: 2.prune.done
	$(ARCH) resolve \
		--in $(CHECKPOINTS)/2.pruned \
		--out $(CHECKPOINTS)/3.resolved \
		--dot 3.resolved.dot \
		--stats 3.resolved.stats \
		> 3.resolve.log
	touch $@





