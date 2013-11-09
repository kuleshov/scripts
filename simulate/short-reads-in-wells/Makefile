NUM_WELLS:=1000
FRAGMENT_LEN:=10000
FRAGMENTS_PER_WELL:=5

REF:=/home/kuleshov/genomes/r_sphaeroides/rs-2.4.corrected.fna

WELLS:=$(shell python -c "for i in xrange($(NUM_WELLS)): print i,")

run:
	make shortreads.fastq

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

$(WELLS:=.shortreads.bfast.fastq): %.shortreads.bfast.fastq: %.fragments.fasta
	dwgsim \
		-e 0.0 \
		-E 0.0 \
		-d 0 \
		-s 50.0 \
		-C 10 \
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
	
shortreads.fastq: $(WELLS:=.shortreads.bfast.fastq)
	cat $^ > $@

shortreads.subsampled.fastq: shortreads.fastq
	awk '{OFS="\t"; getline seq; getline sep; getline qual; print $$0,seq,sep,qual}' $< | shuf | head -n 100000 | awk '{OFS="\n"; print $$1,$$2,$$3,$$4}' > $@.tmp
	mv $@.tmp $@
