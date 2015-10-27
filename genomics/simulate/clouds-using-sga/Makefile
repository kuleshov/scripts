NUM_WELLS:=400
FRAGMENT_LEN:=10000
FRAGMENTS_PER_WELL:=50

REF:=/home/kuleshov/genomes/mock-meta/mock_metagenome.corrected.formated.fasta

WELLS:=$(shell python -c "for i in xrange($(NUM_WELLS)): print i,")

run:
	make reads.assembled

clean:
	-rm *.bed wells.bed.done *.fasta *.fastq *.shortreads.* *.tmp
###############################################################################
## PREPROCESS REFERENCE

renamed.ref: $(REF)
	python bin/rename_contigs.py --ctg $(REF) --out $@.tmp
	mv $@.tmp $@

###############################################################################
## GENERATE FRAGMENTS PER WELL

wells.bed.done: renamed.ref
	python bin/sample_fragments.py --ref $< --len $(FRAGMENT_LEN) --num_fragments $(FRAGMENTS_PER_WELL) --num_wells $(NUM_WELLS)
	touch $@

$(WELLS:=.fragments.bed): %.fragments.bed: wells.bed.done

$(WELLS:=.fragments.fasta): %.fragments.fasta: %.fragments.bed renamed.ref
	bedtools getfasta -fi renamed.ref -bed $< -fo $@.tmp
	mv $@.tmp $@

###############################################################################
## GENERATE SHORT READS IN EACH WELL

$(WELLS:=.shortreads.fasta): %.shortreads.fasta: %.fragments.fasta
	dwgsim \
		-e 0.0 \
		-E 0.0 \
		-d 0 \
		-s 50.0 \
		-C 50 \
		-1 200 \
		-2 0 \
		-r 0.0 \
		-R 0.0 \
		-y 0.0 \
		-c 0 \
		-z 1 \
		-q ] \
		-P well$* \
		$< $*.shortreads

###############################################################################
## ASSEMBLE READS

reads.assembled: $(WELLS:=.shortreads.fasta)
	-mkdir wells
	-rm -rf wells/*
	bash bin/start_assembly.sh $(NUM_WELLS)
	touch $@
