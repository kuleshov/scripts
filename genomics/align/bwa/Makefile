REF:=
READS1:=
READS2:=
BWA:=/home/kuleshov/bin/bwa

all:
	make reads.bam.bai

ref.fasta: $(REF)
	ln -s $< $@

ref.fasta.bwt: ref.fasta
	$(BWA) index $<

reads.bam: $(READS1) $(READS2) ref.fasta.bwt
	$(BWA) mem ref.fasta $(READS1) $(READS2) -t 8 | samtools view -bS - | samtools sort - reads.tmp
	mv reads.tmp.bam $@

reads.bam.bai: reads.bam
	samtools index $<
