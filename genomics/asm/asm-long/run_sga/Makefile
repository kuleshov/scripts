OVERLAP:=2000
SGA_PATH:=/home/kuleshov/bin/asm/sga/bin
READS:=/home/kuleshov/data/mock_meta/even_cov10/reads.formated.fa

clean:
	rm reads.*

run:
	make reads.assembled
	make stats

stats: 
	python bin/measure_contig_stats.py --contigs default-contigs.fa

reads.preprocessed: $(READS)
	$(SGA_PATH)/sga preprocess $< | ~/lib/fastx/fastx_collapser > $@.tmp
	mv $@.tmp $@
	#cp $< $@

reads.preprocessed.bwt: reads.preprocessed
	$(SGA_PATH)/sga index -t 16 -a ropebwt $<

reads.asqg.gz: reads.preprocessed reads.preprocessed.bwt
	$(SGA_PATH)/sga overlap -t 16 -m $(OVERLAP) $<
	# touch $@

reads.assembled: reads.asqg.gz
	$(SGA_PATH)/sga assemble $<
