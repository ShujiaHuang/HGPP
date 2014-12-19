ASM
===

The directory tree:

```
.
|-- input
|   |-- scaffold1.all.hap
|   |-- scaffold1.fa
|   `-- scaffold1.methy-seq.bam
`-- output
    `-- run.sh 

```

1. For the input files in directory 'input/'

The input files for testing ASE are all in this directory

__scaffold1.methy-seq.bam:__ BAM file, created by aligning all the methylation sequencing reads back to the reference genome(YHref)   
__scaffold1.all.hap:__ Phasing result, created by phasing process    
__scaffold1.fa:__ (Reference sequence: YHref) Fasta sequence of scaffold1 in YHref, which assembled by SOAPdenovo2

2. 'output/'

There is a shell script as an example to test the ASM process.
Do it step by step and it will generate all the output files of ASM

