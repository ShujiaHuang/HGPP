# 1.0 Prepare
ls ../input/fosmid_pileup/*.pileup > pileup.list 
perl ../../../bin/phasing/getMarkerLinkageForGenome.pl ../input/scaffold206.wgs.het.snp pileup.list ../input/fosmid_regions.lst > scaffold206.tran.formate
perl ../../../bin/phasing/filterFosmid.pl scaffold206.tran.formate ../input/YH2.depth.abnormal > scaffold206.tran.formate.filtF
perl ../../../bin/phasing/sortByFirstMarker.pl scaffold206.tran.formate.filtF > scaffold206.marker
rm -f scaffold206.tran.formate scaffold206.tran.formate.filtF


# 2.0 Error corection
../../../bin/phasing/ErrCorr -i scaffold206.marker -n 1 -r -o scaffold206.marker.cor 2> log.out
perl ../../../bin/phasing/sortByFirstMarker.pl scaffold206.marker.cor > scaffold206.marker.cor.sort 
# 2.1 Initial connection (hard connection)
../../../bin/phasing/cnn -i scaffold206.marker.cor.sort -n 1 -o scaffold206.marker.cor.sort.cnn 2>> log.out 

rm scaffold206.marker.cor scaffold206.marker.cor.sort # Temporary files, could be delete

# 3.0 Prepare the input files for phasing
mkdir output_readRef output_format_readRef 
perl ../../../bin/phasing/seperate_reads_accrod_overlap.pl scaffold206.marker.cor.sort.cnn ./output_readRef 2 2>> log.out 
perl ../../../bin/phasing/produce_ReFHap_format.pl output_readRef output_format_readRef

# 3.1 RefHap to phase the het-markers
mkdir output_RefHap
java -cp ../../../bin/phasing/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -c 2 -v output_format_readRef/read_1.ref.H output_format_readRef/read_1.ref.F output_RefHap/read_1.ref.phase > output_RefHap/read_1.ref.log 
java -cp ../../../bin/phasing/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -c 2 -v output_format_readRef/read_2.ref.H output_format_readRef/read_2.ref.F output_RefHap/read_2.ref.phase > output_RefHap/read_2.ref.log 
java -cp ../../../bin/phasing/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -c 2 -v output_format_readRef/read_3.ref.H output_format_readRef/read_3.ref.F output_RefHap/read_3.ref.phase > output_RefHap/read_3.ref.log 
java -cp ../../../bin/phasing/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -c 2 -v output_format_readRef/read_4.ref.H output_format_readRef/read_4.ref.F output_RefHap/read_4.ref.phase > output_RefHap/read_4.ref.log 
java -cp ../../../bin/phasing/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -c 2 -v output_format_readRef/read_5.ref.H output_format_readRef/read_5.ref.F output_RefHap/read_5.ref.phase > output_RefHap/read_5.ref.log 
java -cp ../../../bin/phasing/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -c 2 -v output_format_readRef/read_6.ref.H output_format_readRef/read_6.ref.F output_RefHap/read_6.ref.phase > output_RefHap/read_6.ref.log 

# 4.0 Convert the files' format creat by RefHap
mkdir output_convert
perl ../../../bin/phasing/convert_refhap_result.pl output_readRef/ output_format_readRef/ output_RefHap/ output_convert/
# 4.1 Combine all the phasing result
perl ../../../bin/phasing/combine_hap.pl -in output_convert/ -out ./scaffold206.all.hap

# 5.0 Clean 
rm -rf scaffold206.marker.cor.sort.cnn output_readRef output_format_readRef output_RefHap output_convert log.out
