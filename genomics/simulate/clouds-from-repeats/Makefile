NUM_WELLS:=384
FRAGMENT_LEN:=8000
FRAGMENTS_PER_WELL:=500

REF:=/home/kuleshov/data/mock_metagenome.corrected.fasta
REPEATS:=/home/kuleshov/metagenomica/find_repeats_meta/5kb/repeats.grouped.bed

WELLS:=$(shell python -c "for i in xrange($(NUM_WELLS)): print i,")

run:
	make clouds.assembled.bed
	make stats

clean:
	-rm *.bed wells.done *.ids *.tmp

###############################################################################
## PREPROCESS REFERENCE

renamed.ref: $(REF)
	python bin/rename_contigs.py --ctg $(REF) --out $@.tmp
	mv $@.tmp $@

###############################################################################
## GENERATE FRAGMENTS PER WELL

wells.done: renamed.ref
	python bin/sample_fragments.py --ref $< --len $(FRAGMENT_LEN) --num_fragments $(FRAGMENTS_PER_WELL) --num_wells $(NUM_WELLS)
	touch $@

$(WELLS:=.fragments.bed): %.fragments.bed: wells.done

###############################################################################
## SELECT REPEATS THAT OCCUR IN EACH WELL

repeats.all.sorted.bed: $(REPEATS)
	sort -k 4 $< > $@.tmp
	mv $@.tmp $@

$(WELLS:=.repeats.ids): %.repeats.ids: %.fragments.bed repeats.all.sorted.bed
	bedtools intersect -a $< -b repeats.all.sorted.bed -wa -wb | cut -f 7 | sort | uniq -d > $@.tmp
	mv $@.tmp $@

$(WELLS:=.repeats.bed): %.repeats.bed: %.repeats.ids repeats.all.sorted.bed
	join -1 4 -2 1 repeats.all.sorted.bed 0.repeats.ids | cut -d ' ' -f 2- | tr ' ' '\t' > $@.tmp
	mv $@.tmp $@

###############################################################################
## COMPUTE EFFECTIVE CLOUDS

$(WELLS:=.assembled.bed): %.assembled.bed: %.fragments.bed %.repeats.bed 
	bedtools subtract -a $< -b $*.repeats.bed  > $@.tmp
	mv $@.tmp $@

clouds.assembled.bed: $(WELLS:=.assembled.bed)
	cat $^ > $@.tmp
	mv $@.tmp $@

stats:
	python bin/measure_stats.py --contigs clouds.assembled.bed

###############################################################################
## WHICH REPEATS CAN BE RESOLVED?

$(WELLS:=.repeats.counts): %.repeats.counts: %.fragments.bed repeats.all.sorted.bed
	bedtools intersect -a $< -b repeats.all.sorted.bed -wa -wb | cut -f 7 | sort | uniq -c > $@.tmp
	mv $@.tmp $@

repeats.solvability: $(WELLS:=.repeats.counts)
	python count_repeats.py --wells $(NUM_WELLS) --solvability $@.tmp
	mv $@.tmp $@

solvability: repeats.solvability
	cat $< | tail -n +386 | sort -k 3 -r
