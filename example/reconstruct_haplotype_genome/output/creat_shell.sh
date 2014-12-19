prgdir=../../../bin/reconstruct_haplotype_genome
mkdir fastaDirByHap

# extract the FAH sequence covering scaffold206, and seperate them in independent directories. Stored with suffix "*.query.fa"
perl $prgdir/extract_contig_for_genome.pl ../input/scaffold206.list ../input/scaffold206.all.hap ../input/fosmid.fa fastaDirByHap/

############################################################
# Now we're going to creat two shell for the following steps
############################################################

# extract the YHref sequence within scaffold206. Stored with suffix "*.target.fa"
perl $prgdir/extract_target_for_genome.pl ../input/scaffold206.fa ../input/scaffold206.list ../input/scaffold206.all.hap fastaDirByHap/

# generate shell script for step01. As an example, in step01 we only test the hap "scaffold206_15"
perl $prgdir/generate_shell.pl fastaDirByHap/ ./ | grep scaffold206_15 > step01_assemble_for_scaffold206_15.sh

# generate shell script for step02. As an example, in step02 we still only test the hap "scaffold206_15"
mkdir -p RevisedSeq/scaffold206_15
perl $prgdir/get_align_shell.pl fastaDirByHap/scaffold206_15/cor.format.fa ../input/scaffold206_15.pe1.fq.gz ../input/scaffold206_15.pe2.fq.gz RevisedSeq/scaffold206_15 > step02_revise_for_scaffold206_15.sh

