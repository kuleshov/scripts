REPEAT_SIZE:=150
SHIFT:=5
MIN_SCORE:=125
BLAT:=/home/kuleshov/bin/blat
REF:=/home/kuleshov/genomes/r_sphaeroides/rs-2.4.fna

run:
	make measure_stats

measure_stats: contigs.bed
	python bin/measure_stats.py --contigs $<

clean:
	-rm segments.fa out.psl repeat* contigs.bed

renamed.ref: $(REF)
	python bin/rename_contigs.py --ctg $(REF) --out $@.tmp
	mv $@.tmp $@

segments.fa: renamed.ref
	python bin/get_segments.py --len $(REPEAT_SIZE) --shift $(SHIFT) --ref $< --fa $@.tmp
	mv $@.tmp $@

out.psl: renamed.ref segments.fa
	$(BLAT) $^ $@ -minScore=$(MIN_SCORE)

repeat.counts: out.psl
	cat $< | cut -f 10 | uniq -c | sort -n -k 1 -r | grep -v '1 ' > $@

repeats.bed: repeat.counts
	cat $< | sed -e 's/^[ ]*//g' | cut -f 2 -d ' ' | tr '_' '\t' | awk -F '\t' '{ $$2++; $$3++; print $$0}' | tr ' ' '\t' > $@

repeats.merged.bed: repeats.bed
	cat $< | bedtools sort | bedtools merge > $@.tmp
	mv $@.tmp $@

ref.genome: renamed.ref
	python bin/get_genome_size.py --ref $< --genome $@.tmp
	mv $@.tmp $@

contigs.bed: repeats.grouped.bed ref.genome
	bedtools complement -i $< -g ref.genome > $@.tmp
	mv $@.tmp $@

repeats.grouped.bed: out.psl
	python bin/assign_repeat_ids.py --psl $< --repeat_bed $@.tmp
	mv $@.tmp $@
