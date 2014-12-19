# Step 1
ls ../input/bam/*/*.bam > bam.lst
perl ../../../bin/anchor/gen_shell.pl ../../../bin/anchor/anchor_fosmid bam.lst 20 fosmid_regions > step01.find_fosmid.sh

# Step 2 Can only launch after step1 finish
echo "ls $PWD/fosmid_regions/*/*.region > fosmid_regions.lst " > step02.merge_fosmid.sh
echo "perl ../../../bin/anchor/cat.pl fosmid_regions.lst > fosmid.regions.txt" >> step02.merge_fosmid.sh
