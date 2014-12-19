#Author : Shujia Huang
#Date   : 2010 - 04 - 21
#!/usr/bin/perl 

use strict;
use warnings;

my ($file,$index,$aphle) = @ARGV;
my $m = `wc -l $file`;
$m =~ m/(\d+)/;
$m = $1;
my $k = 0;

my ($p, $pre_p);
my $fdr_q;
open (IN,$file) or die "Cannot open : $file\n";
while (<IN>){
	chomp;
	++$k;
	$p = (split (/\t/,$_))[$index-1];
	die "Your input file is not been sorted from samllest to greatest by p_value!\n" if ( $k > 1 && $p < $pre_p ); 

	#$fdr = $k * $aphle / $m;
	$fdr_q = $p * $m / $k;
	if ( $fdr_q < $aphle ){
		print "$_\n";
	}

	$pre_p = $p;
}
close (IN);


