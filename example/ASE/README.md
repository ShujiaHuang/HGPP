ASE
===

The directory tree:

```
.
|-- input   
|   |-- scaffold101.RNA-seq.bam   
|   |-- scaffold101.all.hap       
|   `-- scaffold101.refGene.map2scaffold.gff     
`-- output     
    `-- run.sh 

```

1. For the input files in directory 'input/'

The input files for testing ASE are all in this directory

__scaffold101.refGene.map2scaffold.gff:__ Gene file, created by mapping YHref to hg19 reference   
__scaffold101.RNA-seq.bam:__ BAM file, created by aligning all the RNA-seq reads back to the reference genome(YHref)   
__scaffold101.all.hap:__ Phasing result, created by phasing process    

2. 'output/'

There is a shell script as an example to test the ASE process.
Do it step by step and it will generate all the output files of ASE

