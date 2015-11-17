CTMT:=../spades-shotgun-to-ref/scaffolds.containment
FASTA:=../../asm/spades-shotgun/out/scaffolds.fasta
BAM:=../shotgun-to-spades-shotgun/reads.bam 

overlaps:
	python ./overlaps.py -c $(CTMT) -f $(FASTA)

coarse-pe:
	python ./coarse-pe-connections.py \
	  -b $(BAM) \
	  -f $(FASTA)

explore-pe:
	python ./explore-pe.py \
	  -b $(BAM) \
	  -f $(FASTA)

connections:
	make edges.tsv

edges.tsv:
	python ./pe-connections.py \
	  -b $(BAM) \
	  -f $(FASTA) \
	  -e $@ \
	  -l 125 \
	  -s 200

besst.edges.tsv:
	python ./besst2tsv.py \
	  -s ~/graph.vertices.besst.dump \
	  -e ~/graph.besst.dump \
	  -t besst.edges.tsv
