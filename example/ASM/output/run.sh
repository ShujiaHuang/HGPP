# Step 1. Distinguish the reads from alignment result base on the phasing markers. Use the same program in ASE
#mkdir bam
#../../../extensions/ASE/DistinguishReads -i ../input/scaffold1.methy-seq.bam -f ../input/scaffold1.all.hap -r ../input/scaffold1.fa -B -o bam/scaffold1 && samtools sort bam/scaffold1.1.bam bam/scaffold1.1 && samtools sort bam/scaffold1.2.bam bam/scaffold1.2 && samtools sort bam/scaffold1.homo.bam bam/scaffold1.homo && samtools sort bam/scaffold1.ambiguity.bam bam/scaffold1.ambiguity

# Step 2. Created methylation cout files for the two haplotypes
samtools_path=~/Bin/software_pip/samtools-0.1.19
python ../../../extensions/ASM/methratio.py -c scaffold1 -d ../input/scaffold1.fa -s $samtools_path -o scaffold1.1.cout bam/scaffold1.1.bam && echo "** 1.cout done **" 
python ../../../extensions/ASM/methratio.py -c scaffold1 -d ../input/scaffold1.fa -s $samtools_path -o scaffold1.2.cout bam/scaffold1.2.bam && echo "** 2.cout done **"

# Step 3. Different methylaotion regions(DMR) detection
../../../extensions/ASM/tDMR_detection/tdmr slide -c CG scaffold1.1.cout scaffold1.2.cout > scaffold1.CG && echo "** Window slide done **"
../../../extensions/ASM/tDMR_detection/tdmr extend -b -t 1 -c CG scaffold1.1.cout scaffold1.2.cout scaffold1.CG > scaffold1.CG.ext 

# Step 4. Quality control. FDR
perl ../../../extensions/ASM/unique.pl scaffold1.CG.ext | sort -g -k 4,4 > t.tmp && mv -f t.tmp scaffold1.CG.ext && echo "** extend region done **"
perl ../../../extensions/ASM/filterFDR.pl scaffold1.CG.ext 4 0.01 | awk '$18>=0.9' | sort -k 1,1 -k 2,2n > final.dmr

