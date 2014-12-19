# Step 1. Distinguish the reads from alignment result base on the phasing markers
mkdir bam
../../../extensions/ASE/DistinguishReads -i ../input/scaffold101.RNA-seq.bam -f ../input/scaffold101.all.hap -o bam/scaffold101 && samtools sort bam/scaffold101.1.bam bam/scaffold101.1 && samtools sort bam/scaffold101.2.bam bam/scaffold101.2 && samtools sort bam/scaffold101.homo.bam bam/scaffold101.homo && samtools sort bam/scaffold101.ambiguity.bam bam/scaffold101.ambiguity

# Step 2. Calculting the coverage depth for each gene
ls $PWD/bam/*.1.bam > 1.bam.lst
ls $PWD/bam/*.2.bam > 2.bam.lst
perl ../../../extensions/ASE/calcu_refGene_cover.pl ../input/scaffold101.refGene.map2scaffold.gff scaffold101 1.bam.lst > scaffold101.1.gene.cov 
perl ../../../extensions/ASE/calcu_refGene_cover.pl ../input/scaffold101.refGene.map2scaffold.gff scaffold101 2.bam.lst > scaffold101.2.gene.cov 

# Step 3. Computing the ASE genes
perl ../../../extensions/ASE/ASE.pl scaffold101.1.gene.cov scaffold101.2.gene.cov | awk '$7>0.5 || $7 < -0.5' > ASE.gene

