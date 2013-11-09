REF:=/home/kuleshov/data/mock_metagenome.corrected.fasta
META:=/home/kuleshov/data/mock_metagenome.concentrations
FASTX_PATH:=/home/kuleshov/lib/fastx
LEN:=8000
COV:=10

clean:
	rm reads.*

run:
	make reads.formated.fa
	make reads.formated.qual

###############################################################################
## GENERATE FASTA READS

reads.raw.fa: $(REF)
	bin/sample_random_reads.py \
		--len $(LEN) \
		--cov $(COV) \
		--ref $(REF) \
		--meta $(META) \
		--fasta $@.tmp \
		--qual reads.raw.qual.tmp
	mv $@.tmp $@
	mv reads.raw.qual.tmp reads.raw.qual

reads.formated.fa: reads.raw.fa
	cat $< | $(FASTX_PATH)/fasta_formatter -w 80 > $@.tmp
	mv $@.tmp $@

reads.formated.qual: reads.raw.qual
	cat $< | $(FASTX_PATH)/fasta_formatter -w 80 > $@.tmp
	mv $@.tmp $@

###############################################################################
## SEPARATE READS

reads.cut.fa: reads.raw.fa
	separate_reads.py --raw $< --len 3000 --overlap 2000 > $@.tmp
	mv $@.tmp $@

