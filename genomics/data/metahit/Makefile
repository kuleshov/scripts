N:=30

samples.random.ids: ./MetaHitIDs2EBIacc.txt
	cat $< | shuf | head -n $(N) | cut -f 2 > $@

dl.done: ./samples.random.ids ./samples.txt
	bash dl.sh $< ./samples.txt
	touch $@

count.done: samples.random.ids dl.done
	bash count.sh $<
	touch $@
