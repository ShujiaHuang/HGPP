README
======

1. For directory "input", there are mutiple files, like:  
__scaffold206.list   :__  hap_id and the correspongding fosmid_id, FAH_id, scaffold_id, scaffold_start, scaffold_end
__scaffold206.all.hap:__  The phased results by our phasing pipeline
__fosmid.fa          :__  FAH sequence for scaffold206
__scaffold206.fa     :__  YHref sequence for scaffold206
__scaffold206_15.pe1.fq.gz:__ correspongding reads for hap "scaffold206_15" extracted from the raw reads of fosmid pooling sequencing
__scaffold206_15.pe2.fq.gz:__ correspongding reads for hap "scaffold206_15" extracted from the raw reads of fosmid pooling sequencing

2. For directory "output"

There is a shell script as an example to test the haplotype sequence recontruction process.
Do it step by step and it will generate all the result


