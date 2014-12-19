HGPP
====

1. Overview
-----------

![pipeline plot](/img/pipeline.png)

HGPP(Human Genome haplotype-resolved Program Package) provides a pipeline for building haplotype resolved human genome by fosmid sequencing 

__Contributors:__ Shujia Huang, Yuhui Sun, Xin Tong and Peng Sun <br/>
__Institute   :__ BGI-Shenzhen                                   <br/>
__Last Version:__ 2014-06-18                                     <br/>

This package contains: src, extensions, example and release directories

For `src`: source code for HGPP     
For `img`: Pipeline plot figure    
For `extensions`: Some basic programs for ASE and ASM analysis     
For `example`: Provide some examples for the main processes of HGPP        


2. Install the software
-----------------------

This package is quit easy to install in your directory, just use:

```
make
```

Then all the required programs will all be compiled and installed automatically. You can find the main part executable programs in `bin` and find the programs for ASE/ASM analysis in `extensions`, respectively 

3. How to use HGPP
-------------------

You should prepare all the alignment files([BAM/SAM format](http://samtools.github.io/hts-specs/SAMv1.pdf)) before you get start to run this pipeline.

The processes in HGPP:    

1) Find the fosmid region   
2) Markers error correction     
3) Fosmid initial connection       
4) RefHap phasing process       
5) haplotype genome reconstruction       
6) ASE and ASM       

We have provided some cases for how we use HGPP in our project in directory "example"(Will be created when INSTALL HGPP). You can find the detail tutorial information of HGPP in it. 


4. Please cite the paper
-----------------------

De novo assembly of a haploid-resolved human genome. (Under review)



