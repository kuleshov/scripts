CELPATH:=/home/kuleshov/bin/asm/wgs-7.0/Linux-amd64/bin
# CELPATH:=/home/kuleshov/bin/asm/MaSuRCA-2.0.3/CA/Linux-amd64/bin
INPUTPATH:=/home/kuleshov/data/mock_meta/even_cov10
READS:=$(INPUTPATH)/reads.raw.fa
QUALS:=$(INPUTPATH)/reads.raw.qual

clean:
	-rm reads.fa reads.qual input.frg
	-rm -rf newrun

reads.fa: $(READS)
	# cat $< | head -n 2000 > $@.tmp
	cat $< > $@.tmp
	mv $@.tmp $@

reads.qual: $(QUALS)
	# cat $< | head -n 2000 > $@.tmp
	cat $< > $@.tmp
	mv $@.tmp $@

input.frg: reads.fa reads.qual
	/home/kuleshov/bin/asm/wgs-7.0-merylproblem/Linux-amd64/bin/fastaToCA -l input -s reads.fa -q reads.qual > $@.tmp
	mv $@.tmp $@

run: input.frg
	$(CELPATH)/runCA -d $(CURDIR)/newrun -p test -s run.spec $<
