#Author : Shujia Huang
#Date   : 2011/04/01
#!/usr/bin/perl
use strict;
use warnings;

my ( $sample1_file, $sample2_file, $signif_p ) = @ARGV;
die "perl $0 [sample1] [sample2] > output\n" if ( @ARGV < 2 );

$signif_p ||= 0.01;
my %CALC; #[DO NOT DELETE THIS PARAMETERS] Record the calculated result to avoid calculating again.
my ( %data1, %data2 );
read_file( $sample1_file, \%data1 );
read_file( $sample2_file, \%data2 );

my ( %common, $p_value, $logExp, $upDown, @data );
for my $k ( keys %data1 ) {

	if ( exists $data2{$k} ) {

		next if ( $data1{$k}->[2] < 10 && $data2{$k}->[2] < 10 );

		$p_value = possion_test( $data1{$k}->[2], $data1{$k}->[4], $data2{$k}->[2], $data2{$k}->[4] );

		if ( $data1{$k}->[3] == 0 ) {
			$logExp = "-inf";
		} elsif ( $data2{$k}->[3] == 0 ) {
			$logExp = "+inf";
		} else {
			$logExp = log( $data1{$k}->[3]/$data2{$k}->[3] ) / log(2);
		}

		$upDown = "-"    if ( $data1{$k}->[3] == $data2{$k}->[3] );
		$upDown = "Up"   if ( $data1{$k}->[3]  < $data2{$k}->[3] );
		$upDown = "Down" if ( $data1{$k}->[3]  > $data2{$k}->[3] );
		# scaffold100 BAG4 90 12.9560355994579 4940661
		push ( @data, [ $data1{$k}->[0],$data1{$k}->[1],$data1{$k}->[2],$data2{$k}->[2],$data1{$k}->[3],$data2{$k}->[3],$logExp,$upDown,$p_value ] ) if ( $p_value < $signif_p );

		$common{$k} = 1;
	}
}

my $m = scalar @data;
@data = sort { $a->[8] <=> $b->[8] } @data;

print "#Scaffold\tGene_Symbol\tExp1\tExp2\tRPKM1\tRPKM2\tlog2Ratio(RPKM1/RPKM2)\tUp-Down-Regular\tP-value\tFDR\n";
for ( my $i = 0; $i < $m; ++$i ) {

	my $fdr_q = $data[$i]->[8] * $m / ( $i + 1 );
	if ( $fdr_q < $signif_p ) {
		print join "\t", ( @{ $data[$i] }, $fdr_q ), "\n";
	}
}


###############################
sub read_file {
	my ( $file, $data ) = @_;

	my ( $start, $end );
	open ( I, $file ) or die "Cannot open file : $file\n";
	while (<I>){ #
		
		next if /^#/;
		chomp;
		#scaffold100     BAG4    1406    90      12.9560355994579        4940661
		my ( $contig, $gene, $tag1, $exp, $total_tag ) = ((split)[0,1,3..5]); # scaffold100 BAG4 90 12.9560355994579 4940661
		my $key      = $contig.":".$gene;
		$$data{$key} = [ $contig, $gene, $tag1, $exp, $total_tag ];
	}
	close ( I );
}

sub possion_test{

	my ( $x,$n1, $y,$n2 ) = @_;

	#die "[ERROR]The number should not be 0 in possion_test() !\n" if ( $n1 == 0 || $n2 == 0 );

	my ($log_p_value, $p_value) = ( 0, 0 );
	if ( $n1 ) {
		$log_p_value = $y * ( log($n2) - log($n1) ) + factorial_loge( $x+$y ) - 
					   ( factorial_loge($x) + factorial_loge($y) + ($x+$y+1) * log( 1 + $n2 / $n1) );
	    $p_value    += exp( $log_p_value ); # The probability of the gene's expression is equally in these two samples.
	} else {
		$p_value     = 0;
	}

	return $p_value;
}

sub factorial_loge {
	my ( $num ) = @_;
	my $result  = 0;

	return $CALC{$num} if ( defined $CALC{$num} );

	if ( $num > 1 ){

		for ( my $i = 1; $i <= $num; ++$i ){

			$result += log($i);
		}	
	} else {

		if ( $num < 0 ) {
			die "[ERROR]Negative data when calculte factorial\n";
		} else {
			$result = 0;
		}
	}
	$CALC{$num} = $result;
	return $result;
}

sub swap {
	my ( $dat1, $dat2 ) = @_;
	my $tmp = $$dat2;
	$$dat2  = $$dat1;
	$$dat1  = $tmp;
}








